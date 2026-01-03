package OperatorAlgebra

import (
	"GoProject/gridData"
	"errors"
	"fmt"
	"math"
	"math/cmplx"

	"gonum.org/v1/gonum/mat"
)

// KeDvrBasis represents the kinetic energy operator
type KeDvrBasis struct {
	grid       *gridData.RadGrid
	mass       float64
	ndims      int
	dx2        float64
	invMassDx2 float64
	diagTerm   float64
	kMat       *mat.Dense

	kMatCached bool
	eigCached  bool
	eigVals    []float64
	eigVecs    *mat.Dense
}

// NewKeDVR creates and initializes DVR basis
func NewKeDVR(grid *gridData.RadGrid, mass float64) *KeDvrBasis {
	ndim := int(grid.NPoints())
	dx2 := grid.DeltaR() * grid.DeltaR()
	invMassDx2 := 1.0 / (mass * dx2)
	diagTerm := math.Pi * math.Pi / (6.0 * mass * dx2)

	return &KeDvrBasis{
		grid:       grid,
		mass:       mass,
		ndims:      ndim,
		dx2:        dx2,
		invMassDx2: invMassDx2,
		diagTerm:   diagTerm,
		kMat:       mat.NewDense(ndim, ndim, nil),
		kMatCached: false,
		eigCached:  false,
		eigVals:    nil,
		eigVecs:    nil,
	}
}

// isZeroToInfinity determines if the grid domain is [0, infinity)
func (k *KeDvrBasis) isZeroToInfinity() bool {
	return math.Abs(k.grid.RMax()+k.grid.RMin()) >= 1
}

// MinInftyToInfty constructs the kinetic energy matrix
func (k *KeDvrBasis) MinInftyToInfty(m *mat.Dense) {
	for i := 0; i < k.ndims; i++ {
		m.Set(i, i, k.diagTerm)
	}
	for i := 1; i < k.ndims; i++ {
		for j := 0; j < i; j++ {
			diff := i - j
			sign := float64(1 - 2*(diff&1))
			val := sign * k.invMassDx2 / float64(diff*diff)
			m.Set(i, j, val)
			m.Set(j, i, val)
		}
	}
}

// ZeroToInfinity constructs the kinetic energy matrix for [0, ∞) domain
func (k *KeDvrBasis) ZeroToInfinity(m *mat.Dense) {
	constVal := 0.25 * k.invMassDx2
	for i := 0; i < k.ndims; i++ {
		m.Set(i, i, k.diagTerm-constVal/float64((i+1)*(i+1)))
	}

	for i := 1; i < k.ndims; i++ {
		for j := 0; j < i; j++ {
			diff := i - j
			add := i + j
			sign := float64(1 - 2*(diff&1))
			diffSq := float64(diff * diff)
			addSq := float64(add * add)
			val := sign * k.invMassDx2 * (1.0/diffSq - 1.0/addSq)
			m.Set(i, j, val)
			m.Set(j, i, val)
		}
	}
}

// GetMat returns the kinetic energy matrix, using cache if available
func (k *KeDvrBasis) GetMat() *mat.Dense {
	if !k.kMatCached {
		if k.isZeroToInfinity() {
			k.ZeroToInfinity(k.kMat)
		} else {
			k.MinInftyToInfty(k.kMat)
		}
		k.kMatCached = true
		k.eigCached = false
	}
	return k.kMat
}

func (k *KeDvrBasis) ensureEigen() error {
	if k.eigCached {
		return nil
	}

	KM := k.GetMat()

	// Convert Dense to Symmetric
	sym := mat.NewSymDense(k.ndims, nil)
	for i := 0; i < k.ndims; i++ {
		for j := i; j < k.ndims; j++ {
			sym.SetSym(i, j, KM.At(i, j))
		}
	}

	var es mat.EigenSym
	ok := es.Factorize(sym, true)
	if !ok {
		return errors.New("eigendecomposition failed (EigenSym.Factorize returned false)")
	}

	k.eigVals = es.Values(nil)
	k.eigVecs = mat.NewDense(k.ndims, k.ndims, nil)
	es.VectorsTo(k.eigVecs)

	k.eigCached = true
	return nil
}

// ExpDtTo applies exp(Dt * K) to In and writes to Out (real -> real)
func (k *KeDvrBasis) ExpDtTo(Dt float64, In []float64, Out []float64) error {
	if len(In) != k.ndims || len(Out) != k.ndims {
		return fmt.Errorf("vector length mismatch: got In=%d Out=%d, expected %d", len(In), len(Out), k.ndims)
	}
	if err := k.ensureEigen(); err != nil {
		return err
	}

	V := k.eigVecs
	lams := k.eigVals

	// tmp = V^T * In  (real)
	tmp := make([]float64, k.ndims)
	for i := 0; i < k.ndims; i++ {
		sum := 0.0
		for j := 0; j < k.ndims; j++ {
			sum += V.At(j, i) * In[j] // V^T: (i,j) = V(j,i)
		}
		tmp[i] = sum
	}

	// tmp[i] *= exp(Dt * lambda_i)
	for i := 0; i < k.ndims; i++ {
		tmp[i] *= math.Exp(Dt * lams[i])
	}

	// Out = V * tmp
	for i := 0; i < k.ndims; i++ {
		sum := 0.0
		for j := 0; j < k.ndims; j++ {
			sum += V.At(i, j) * tmp[j]
		}
		Out[i] = sum
	}

	return nil
}

// ExpIdtTo applies exp(i * Dt * K) to In and writes to Out (complex -> complex)
func (k *KeDvrBasis) ExpIdtTo(Dt float64, In []complex128, Out []complex128) error {
	if len(In) != k.ndims || len(Out) != k.ndims {
		return fmt.Errorf("vector length mismatch: got In=%d Out=%d, expected %d", len(In), len(Out), k.ndims)
	}
	if err := k.ensureEigen(); err != nil {
		return err
	}

	V := k.eigVecs
	lams := k.eigVals

	// tmp = V^T * In  (complex)
	tmp := make([]complex128, k.ndims)
	for i := 0; i < k.ndims; i++ {
		var sum complex128
		for j := 0; j < k.ndims; j++ {
			sum += complex(V.At(j, i), 0) * In[j]
		}
		tmp[i] = sum
	}

	// tmp[i] *= exp(i * Dt * lambda_i)
	for i := 0; i < k.ndims; i++ {
		tmp[i] *= cmplx.Exp(complex(0, Dt*lams[i]))
	}

	// Out = V * tmp
	for i := 0; i < k.ndims; i++ {
		var sum complex128
		for j := 0; j < k.ndims; j++ {
			sum += complex(V.At(i, j), 0) * tmp[j]
		}
		Out[i] = sum
	}

	return nil
}

// ExpDt computes exp(dt * K) and applies it to input vector (allocating output)
func (k *KeDvrBasis) ExpDt(dt float64, in []float64) ([]float64, error) {
	if len(in) != k.ndims {
		return nil, fmt.Errorf("input vector length %d doesn't match basis dimension %d", len(in), k.ndims)
	}
	out := make([]float64, k.ndims)
	if err := k.ExpDtTo(dt, in, out); err != nil {
		return nil, err
	}
	return out, nil
}

// ExpDtInPlace computes exp(dt * K) and applies it to the input vector in-place
func (k *KeDvrBasis) ExpDtInPlace(dt float64, inOut []float64) error {
	if len(inOut) != k.ndims {
		return fmt.Errorf("vector length %d doesn't match basis dimension %d", len(inOut), k.ndims)
	}
	tmp := make([]float64, k.ndims)
	if err := k.ExpDtTo(dt, inOut, tmp); err != nil {
		return err
	}
	copy(inOut, tmp)
	return nil
}

// ExpIdt computes exp(i*dt*K) and applies it to a complex input vector (allocating output)
func (k *KeDvrBasis) ExpIdt(dt float64, in []complex128) ([]complex128, error) {
	if len(in) != k.ndims {
		return nil, fmt.Errorf("input vector length %d doesn't match basis dimension %d", len(in), k.ndims)
	}
	out := make([]complex128, k.ndims)
	if err := k.ExpIdtTo(dt, in, out); err != nil {
		return nil, err
	}
	return out, nil
}

func (k *KeDvrBasis) ExpIdtInPlace(dt float64, inOut []complex128) error {
	if len(inOut) != k.ndims {
		return fmt.Errorf("vector length %d doesn't match basis dimension %d", len(inOut), k.ndims)
	}
	tmp := make([]complex128, k.ndims)
	if err := k.ExpIdtTo(dt, inOut, tmp); err != nil {
		return err
	}
	copy(inOut, tmp)
	return nil
}

// RealDiagonalize diagonalizes the real kinetic energy matrix (fixed)
func (k *KeDvrBasis) RealDiagonalize() (eigenvalues []float64, eigenvectors *mat.Dense, err error) {
	if err = k.ensureEigen(); err != nil {
		return nil, nil, err
	}

	eigenvalues = make([]float64, k.ndims)
	copy(eigenvalues, k.eigVals)

	eigenvectors = mat.NewDense(k.ndims, k.ndims, nil)
	eigenvectors.Copy(k.eigVecs)
	return eigenvalues, eigenvectors, nil
}

func (k *KeDvrBasis) Clone() *KeDvrBasis {
	newK := &KeDvrBasis{
		grid:       k.grid,
		mass:       k.mass,
		ndims:      k.ndims,
		dx2:        k.dx2,
		invMassDx2: k.invMassDx2,
		diagTerm:   k.diagTerm,
		kMat:       mat.NewDense(k.ndims, k.ndims, nil),
		kMatCached: false,

		eigCached: false,
		eigVals:   nil,
		eigVecs:   nil,
	}

	if k.kMatCached {
		newK.kMat.Copy(k.kMat)
		newK.kMatCached = true
	}

	// If eigen cache exists, copy it too (optional but handy)
	if k.eigCached && k.eigVals != nil && k.eigVecs != nil {
		newK.eigVals = make([]float64, len(k.eigVals))
		copy(newK.eigVals, k.eigVals)
		newK.eigVecs = mat.NewDense(k.ndims, k.ndims, nil)
		newK.eigVecs.Copy(k.eigVecs)
		newK.eigCached = true
	}

	return newK
}

func (k *KeDvrBasis) Ndims() int {
	return k.ndims
}

func (k *KeDvrBasis) Mass() float64 {
	return k.mass
}

func (k *KeDvrBasis) Clear() {
	k.kMatCached = false
	k.eigCached = false
}
