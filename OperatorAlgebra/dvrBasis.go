package OperatorAlgebra

import (
	"GoProject/gridData"
	"errors"
	"fmt"
	"math"

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
	}
}

func (k *KeDvrBasis) String() string {
	return fmt.Sprintf("KeDVR(xmin:\t%14.7f\n, "+
		"xmax:\t %14.7f\n, "+
		"mass:\t %14.7f\n,"+
		" Dim:%d\n)", k.grid.RMin(), k.grid.RMax(), k.mass, k.ndims)
}

func (k *KeDvrBasis) Ndims() int {
	return k.ndims
}

func (k *KeDvrBasis) Mass() float64 {
	return k.mass
}

func (k *KeDvrBasis) Clear() {
	k.kMatCached = false
}

func (k *KeDvrBasis) isZeroToInfinity() bool {
	return math.Abs(k.grid.RMax()+k.grid.RMin()) >= 1
}

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

func (k *KeDvrBasis) SetMat() {
	if !k.kMatCached {
		if k.isZeroToInfinity() {
			k.ZeroToInfinity(k.kMat)
		} else {
			k.MinInftyToInfty(k.kMat)
		}
		k.kMatCached = true
	}
}

func (k *KeDvrBasis) GetMat() *mat.Dense {
	k.SetMat()
	return k.kMat
}

func (k *KeDvrBasis) ensureEigen(eigvals []float64, eigvecs *mat.Dense) error {
	if eigvecs == nil || eigvals == nil {
		eigvals = make([]float64, k.ndims)
		eigvecs = mat.NewDense(k.ndims, k.ndims, nil)
	}

	if len(eigvals) != k.ndims || eigvecs.Dims() != k.GetMat().Dims() {
		return errors.New("eigenvalues array length mismatch")
	}

	KM := k.GetMat()
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

	eigvals = es.Values(nil)
	eigvecs = mat.NewDense(k.ndims, k.ndims, nil)
	es.VectorsTo(eigvecs)

	return nil
}

func (k *KeDvrBasis) RealDiagonalize(eigenvalues []float64, eigenvectors *mat.Dense) error {
	err := k.ensureEigen(eigenvalues, eigenvectors)
	if err != nil {
		return err
	}
	return nil
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
	}

	if k.kMatCached {
		newK.kMat.Copy(k.kMat)
		newK.kMatCached = true
	}
	return newK
}

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
