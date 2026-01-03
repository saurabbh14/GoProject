package OperatorAlgebra

import (
	"GoProject/gridData"
	"fmt"
	"log"
	"math"
	"math/rand/v2"
	"runtime"

	"github.com/jvlmdr/go-fftw/fftw"
	"gonum.org/v1/gonum/mat"
)

// FourierBasis represents the kinetic energy operator in Fourier/DVR basis
type FourierBasis struct {
	grid     *gridData.RadGrid
	mass     float64
	nPoints  int
	fftPlan  fftw.Plan
	ifftPlan fftw.Plan
	Buff     *fftw.Array
	kValues  []float64
	keValues []float64
}

// FFTInit creates a new Fourier struct for fast fourier transform
func FFTInit(grid *gridData.RadGrid, mass float64) *FourierBasis {
	fb := &FourierBasis{}
	fb.initialize(grid, mass)
	runtime.SetFinalizer(fb, (*FourierBasis).destroy)
	return fb
}

func (f *FourierBasis) Redefine(grid *gridData.RadGrid, mass float64) {
	(*fftw.Plan).Destroy(&f.fftPlan)
	(*fftw.Plan).Destroy(&f.ifftPlan)
	f.initialize(grid, mass)
}

func (f *FourierBasis) initialize(grid *gridData.RadGrid, mass float64) {
	gridPoints := int(grid.NPoints())
	kVal := grid.KValues()

	kValues := make([]float64, gridPoints)
	keValues := make([]float64, gridPoints)
	invTwoMass := 1.0 / (2.0 * mass)

	for i := 0; i < gridPoints; i++ {
		k := kVal[i]
		kValues[i] = k
		keValues[i] = k * k * invTwoMass
	}

	buff := fftw.NewArray(gridPoints)

	for i := 0; i < gridPoints; i++ {
		Re := rand.NormFloat64()
		Im := rand.NormFloat64()
		invMag := 1.0 / math.Sqrt(Re*Re+Im*Im)
		buff.Set(i, complex(Re*invMag, Im*invMag))
	}

	fftPlan := *fftw.NewPlan(buff, buff, fftw.Forward, fftw.Estimate)
	ifftPlan := *fftw.NewPlan(buff, buff, fftw.Backward, fftw.Estimate)

	f.grid = grid
	f.mass = mass
	f.nPoints = gridPoints
	f.fftPlan = fftPlan
	f.ifftPlan = ifftPlan
	f.Buff = buff
	f.kValues = kValues
	f.keValues = keValues
}

func (f *FourierBasis) forwardBuff(in []complex128) {
	copy(f.Buff.Elems, in)
	f.fftPlan.Execute()
}

func (f *FourierBasis) forward(in []complex128, out []complex128) {
	copy(f.Buff.Elems, in)
	f.fftPlan.Execute()
	copy(out, f.Buff.Elems)
}

func (f *FourierBasis) forwardInPlace(InOut []complex128) {
	copy(f.Buff.Elems, InOut)
	f.fftPlan.Execute()
	copy(InOut, f.Buff.Elems)
}

func (f *FourierBasis) backwardBuff(in []complex128) {
	copy(f.Buff.Elems, in)
	f.ifftPlan.Execute()
}

func (f *FourierBasis) backward(in []complex128, out []complex128) {
	copy(f.Buff.Elems, in)
	f.ifftPlan.Execute()
	copy(out, f.Buff.Elems)
}

func (f *FourierBasis) backwardInPlace(InOut []complex128) {
	copy(f.Buff.Elems, InOut)
	f.ifftPlan.Execute()
	copy(InOut, f.Buff.Elems)
}

func (f *FourierBasis) operatorOp(InOut []complex128, Op []float64) {
	if len(f.Buff.Elems) != len(Op) {
		panic(fmt.Sprintf("length mismatch: Buff.Elems=%d, kValues=%d",
			len(f.Buff.Elems), len(Op)))
	}

	copy(f.Buff.Elems, InOut)
	f.fftPlan.Execute()

	for i := range f.Buff.Elems {
		f.Buff.Elems[i] *= complex(Op[i], 0)
	}

	f.ifftPlan.Execute()
	copy(InOut, f.Buff.Elems)
}

func (f *FourierBasis) MomentumOpInPlace(InOut []complex128) {
	f.operatorOp(InOut, f.kValues)
}

func (f *FourierBasis) MomentumOp(In []complex128, Out []complex128) {
	copy(Out, In)
	f.MomentumOpInPlace(Out)
	copy(Out, f.Buff.Elems)
}

func (f *FourierBasis) LaplacianOpInPlace(InOut []complex128) {
	f.operatorOp(InOut, f.keValues)
}

func (f *FourierBasis) LaplacianOp(In []complex128, Out []complex128) {
	copy(Out, In)
	f.LaplacianOpInPlace(Out)
	copy(Out, f.Buff.Elems)
}

func (f *FourierBasis) BuildHamiltonianFGH() *mat.SymDense {
	Tr := make([]float64, f.nPoints)
	pref := 1. / (2.0 * f.mass)
	for r := 1; r <= f.nPoints; r++ {
		k := 2.0 * math.Pi * float64(r) / f.grid.Length()
		Tr[r] = pref * k * k
	}

	H := mat.NewSymDense(f.nPoints, nil)
	for i := 0; i < f.nPoints; i++ {
		for j := 0; j <= i; j++ {
			d := i - j
			sum := 0.0
			for r := 1; r <= f.nPoints; r++ {
				angle := 2.0 * math.Pi * float64(r*d) / float64(f.nPoints)
				sum += math.Cos(angle) * Tr[r]
			}
			hij := (2.0 / float64(f.nPoints)) * sum
			H.SetSym(i, j, hij)
		}
	}
	return H
}

func (f *FourierBasis) DiagonalizeSym() (vals []float64, vecs *mat.Dense) {
	H := f.BuildHamiltonianFGH()
	var es mat.EigenSym
	ok := es.Factorize(H, true)
	if !ok {
		log.Fatal("eigensolve failed")
	}
	vals = es.Values(nil)
	vecs = mat.NewDense(len(vals), len(vals), nil)
	es.VectorsTo(vecs)
	return vals, vecs
}

func (f *FourierBasis) destroy() {
	if f != nil {
		(*fftw.Plan).Destroy(&f.fftPlan)
		(*fftw.Plan).Destroy(&f.ifftPlan)
	}
	fmt.Println("Destroyed")
}

// Clean cleans up FFTW
func (f *FourierBasis) Clean() {
	runtime.SetFinalizer(f, (*FourierBasis).destroy)
	f.destroy()
	fmt.Println("Cleaned")
}
