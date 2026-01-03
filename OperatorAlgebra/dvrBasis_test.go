package OperatorAlgebra

import (
	"GoProject/gridData"
	"math"
	"math/cmplx"
	"testing"

	"gonum.org/v1/gonum/mat"
)

const tolerance = 1e-10

// Helper function to check if two floats are approximately equal
func floatEqual(a, b, tol float64) bool {
	return math.Abs(a-b) < tol
}

// Helper function to check if two complex numbers are approximately equal
func complexEqual(a, b complex128, tol float64) bool {
	return cmplx.Abs(a-b) < tol
}

func TestNewKeDVR(t *testing.T) {
	tests := []struct {
		name    string
		rMin    float64
		rMax    float64
		nPoints int
		mass    float64
	}{
		{"small_symmetric_grid", -5.0, 5.0, 10, 1.0},
		{"zero_to_infinity_grid", 0.0, 10.0, 20, 1.0},
		{"heavy_mass", -10.0, 10.0, 15, 10.0},
		{"light_mass", 0.0, 5.0, 8, 0.5},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			grid, _ := gridData.NewRGrid(tt.rMin, tt.rMax, uint32(tt.nPoints))
			ke := NewKeDVR(grid, tt.mass)

			if ke == nil {
				t.Fatal("NewKeDVR returned nil")
			}

			if ke.Ndims() != tt.nPoints {
				t.Errorf("Ndims() = %d, want %d", ke.Ndims(), tt.nPoints)
			}

			if ke.Mass() != tt.mass {
				t.Errorf("Mass() = %f, want %f", ke.Mass(), tt.mass)
			}
		})
	}
}

func TestIsZeroToInfinity(t *testing.T) {
	tests := []struct {
		name     string
		rMin     float64
		rMax     float64
		expected bool
	}{
		{"symmetric_grid", -5.0, 5.0, false},
		{"zero_to_positive", 0.0, 10.0, true},
		{"small_asymmetry", -0.4, 0.5, false},
		{"large_positive", 1.0, 10.0, true},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			grid, _ := gridData.NewRGrid(tt.rMin, tt.rMax, 10)
			ke := NewKeDVR(grid, 1.0)

			result := ke.isZeroToInfinity()
			if result != tt.expected {
				t.Errorf("isZeroToInfinity() = %v, want %v (rMin=%f, rMax=%f, sum=%f)",
					result, tt.expected, tt.rMin, tt.rMax, tt.rMin+tt.rMax)
			}
		})
	}
}

func TestGetMatSymmetry(t *testing.T) {
	tests := []struct {
		name    string
		rMin    float64
		rMax    float64
		nPoints int
	}{
		{"symmetric_domain", -5.0, 5.0, 10},
		{"zero_to_infinity", 0.0, 10.0, 15},
		{"large_grid", -20.0, 20.0, 50},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			grid, _ := gridData.NewRGrid(tt.rMin, tt.rMax, uint32(tt.nPoints))
			ke := NewKeDVR(grid, 1.0)
			kMat := ke.GetMat()

			if !isSymmetric(kMat) {
				t.Error("Kinetic energy matrix is not symmetric")
			}
		})
	}
}

func TestGetMatCaching(t *testing.T) {
	grid, _ := gridData.NewRGrid(-5.0, 5.0, 10)
	ke := NewKeDVR(grid, 1.0)

	// The first call should compute the matrix
	mat1 := ke.GetMat()

	// Second call should return cached matrix
	mat2 := ke.GetMat()

	// Should be the same pointer
	if mat1 != mat2 {
		t.Error("GetMat() should return cached matrix on second call")
	}
}

func TestClear(t *testing.T) {
	grid, _ := gridData.NewRGrid(-5.0, 5.0, 10)
	ke := NewKeDVR(grid, 1.0)

	// Build cache
	ke.GetMat()
	err := ke.ensureEigen()
	if err != nil {
		return
	}

	// Clear cache
	ke.Clear()

	// Verify the internal state is cleared (indirectly by checking it rebuilds)
	mat1 := ke.GetMat()
	if mat1 == nil {
		t.Error("GetMat() should work after Clear()")
	}
	return
}

func TestMinInftyToInfty(t *testing.T) {
	grid, _ := gridData.NewRGrid(-5.0, 5.0, 5)
	ke := NewKeDVR(grid, 1.0)

	m := mat.NewDense(5, 5, nil)
	ke.MinInftyToInfty(m)

	// Check diagonal is positive (kinetic energy is positive definite)
	for i := 0; i < 5; i++ {
		if m.At(i, i) <= 0 {
			t.Errorf("Diagonal element [%d,%d] = %f should be positive", i, i, m.At(i, i))
		}
	}

	// Check symmetry
	if !isSymmetric(m) {
		t.Error("MinInftyToInfty matrix is not symmetric")
	}

	// Check off-diagonal alternating signs
	for i := 1; i < 5; i++ {
		diff := i - 0
		expectedSign := float64(1 - 2*(diff&1))
		actualSign := math.Copysign(1.0, m.At(i, 0))
		if expectedSign != actualSign {
			t.Errorf("Off-diagonal sign at [%d,0] incorrect: expected %f, got %f",
				i, expectedSign, actualSign)
		}
	}
}

func TestZeroToInfinity(t *testing.T) {
	grid, _ := gridData.NewRGrid(0.0, 10.0, 5)
	ke := NewKeDVR(grid, 1.0)

	m := mat.NewDense(5, 5, nil)
	ke.ZeroToInfinity(m)

	if !isSymmetric(m) {
		t.Error("ZeroToInfinity matrix is not symmetric")
	}

	for i := 0; i < 5; i++ {
		if m.At(i, i) <= 0 {
			t.Errorf("Diagonal element [%d,%d] = %f should be positive", i, i, m.At(i, i))
		}
	}
}

func TestRealDiagonalize(t *testing.T) {
	grid, _ := gridData.NewRGrid(-5.0, 5.0, 10)
	ke := NewKeDVR(grid, 1.0)

	eigenvalues, eigenvectors, err := ke.RealDiagonalize()
	if err != nil {
		t.Fatalf("RealDiagonalize failed: %v", err)
	}

	if len(eigenvalues) != 10 {
		t.Errorf("Expected 10 eigenvalues, got %d", len(eigenvalues))
	}

	r, c := eigenvectors.Dims()
	if r != 10 || c != 10 {
		t.Errorf("Eigenvector matrix dimensions = (%d,%d), want (10,10)", r, c)
	}

	// Eigenvalues should be non-negative for kinetic energy
	for i, ev := range eigenvalues {
		if ev < -tolerance {
			t.Errorf("Eigenvalue[%d] = %f is negative", i, ev)
		}
	}

	// Verify V^T * K * V ≈ diagonal(eigenvalues)
	K := ke.GetMat()
	V := eigenvectors

	// Compute K * V
	KV := mat.NewDense(10, 10, nil)
	KV.Mul(K, V)

	// Compute V^T * K * V
	VtKV := mat.NewDense(10, 10, nil)
	VtKV.Mul(V.T(), KV)

	// Check that it's approximately diagonal with eigenvalues
	for i := 0; i < 10; i++ {
		for j := 0; j < 10; j++ {
			expected := 0.0
			if i == j {
				expected = eigenvalues[i]
			}
			if !floatEqual(VtKV.At(i, j), expected, 1e-8) {
				t.Errorf("V^T*K*V[%d,%d] = %f, expected %f", i, j, VtKV.At(i, j), expected)
			}
		}
	}
}

func TestExpDtTo(t *testing.T) {
	grid, _ := gridData.NewRGrid(-5.0, 5.0, 5)
	ke := NewKeDVR(grid, 1.0)

	in := []float64{1.0, 0.0, 0.0, 0.0, 0.0}
	out := make([]float64, 5)

	err := ke.ExpDtTo(0.0, in, out)
	if err != nil {
		t.Fatalf("ExpDtTo failed: %v", err)
	}

	// With dt=0, exp(0*K) = I, so out should equal in
	for i := range in {
		if !floatEqual(out[i], in[i], tolerance) {
			t.Errorf("ExpDtTo with dt=0: out[%d] = %f, want %f", i, out[i], in[i])
		}
	}
}

func TestExpDtToVectorLengthMismatch(t *testing.T) {
	grid, _ := gridData.NewRGrid(-5.0, 5.0, 5)
	ke := NewKeDVR(grid, 1.0)

	in := []float64{1.0, 0.0, 0.0} // Wrong size
	out := make([]float64, 5)

	err := ke.ExpDtTo(0.1, in, out)
	if err == nil {
		t.Error("ExpDtTo should return error for mismatched vector length")
	}
}

func TestExpIdtTo(t *testing.T) {
	grid, _ := gridData.NewRGrid(-5.0, 5.0, 5)
	ke := NewKeDVR(grid, 1.0)

	in := []complex128{1.0, 0.0, 0.0, 0.0, 0.0}
	out := make([]complex128, 5)

	err := ke.ExpIdtTo(0.0, in, out)
	if err != nil {
		t.Fatalf("ExpIdtTo failed: %v", err)
	}

	// With dt=0, exp(i*0*K) = I, so out should equal in
	for i := range in {
		if !complexEqual(out[i], in[i], tolerance) {
			t.Errorf("ExpIdtTo with dt=0: out[%d] = %v, want %v", i, out[i], in[i])
		}
	}
}

func TestExpIdtUnitarity(t *testing.T) {
	grid, _ := gridData.NewRGrid(-5.0, 5.0, 8)
	ke := NewKeDVR(grid, 1.0)

	// Create a normalized input vector
	in := make([]complex128, 8)
	norm := 0.0
	for i := range in {
		in[i] = complex(float64(i+1), float64(i))
		norm += real(in[i])*real(in[i]) + imag(in[i])*imag(in[i])
	}
	norm = math.Sqrt(norm)
	for i := range in {
		in[i] /= complex(norm, 0)
	}

	out := make([]complex128, 8)
	err := ke.ExpIdtTo(0.1, in, out)
	if err != nil {
		t.Fatalf("ExpIdtTo failed: %v", err)
	}

	// Check that output is also normalized (unitary evolution preserves norm)
	outNorm := 0.0
	for _, v := range out {
		outNorm += real(v)*real(v) + imag(v)*imag(v)
	}
	outNorm = math.Sqrt(outNorm)

	if !floatEqual(outNorm, 1.0, 1e-8) {
		t.Errorf("Unitary evolution did not preserve norm: got %f, want 1.0", outNorm)
	}
}

func TestExpDt(t *testing.T) {
	grid, _ := gridData.NewRGrid(-5.0, 5.0, 5)
	ke := NewKeDVR(grid, 1.0)

	in := []float64{1.0, 0.5, 0.0, -0.5, -1.0}
	out, err := ke.ExpDt(0.01, in)
	if err != nil {
		t.Fatalf("ExpDt failed: %v", err)
	}

	if len(out) != 5 {
		t.Errorf("ExpDt output length = %d, want 5", len(out))
	}
}

func TestExpDtInPlace(t *testing.T) {
	grid, _ := gridData.NewRGrid(-5.0, 5.0, 5)
	ke := NewKeDVR(grid, 1.0)

	in := []float64{1.0, 0.5, 0.0, -0.5, -1.0}
	expected, _ := ke.ExpDt(0.01, in)

	inOut := []float64{1.0, 0.5, 0.0, -0.5, -1.0}
	err := ke.ExpDtInPlace(0.01, inOut)
	if err != nil {
		t.Fatalf("ExpDtInPlace failed: %v", err)
	}

	for i := range expected {
		if !floatEqual(inOut[i], expected[i], tolerance) {
			t.Errorf("ExpDtInPlace[%d] = %f, want %f", i, inOut[i], expected[i])
		}
	}
}

func TestExpIdt(t *testing.T) {
	grid, _ := gridData.NewRGrid(-5.0, 5.0, 5)
	ke := NewKeDVR(grid, 1.0)

	in := []complex128{1.0, 0.5i, 0.0, -0.5i, -1.0}
	out, err := ke.ExpIdt(0.01, in)
	if err != nil {
		t.Fatalf("ExpIdt failed: %v", err)
	}

	if len(out) != 5 {
		t.Errorf("ExpIdt output length = %d, want 5", len(out))
	}
}

func TestExpIdtInPlace(t *testing.T) {
	grid, _ := gridData.NewRGrid(-5.0, 5.0, 5)
	ke := NewKeDVR(grid, 1.0)

	in := []complex128{1.0, 0.5i, 0.0, -0.5i, -1.0}
	expected, _ := ke.ExpIdt(0.01, in)

	inOut := []complex128{1.0, 0.5i, 0.0, -0.5i, -1.0}
	err := ke.ExpIdtInPlace(0.01, inOut)
	if err != nil {
		t.Fatalf("ExpIdtInPlace failed: %v", err)
	}

	// Results should match
	for i := range expected {
		if !complexEqual(inOut[i], expected[i], tolerance) {
			t.Errorf("ExpIdtInPlace[%d] = %v, want %v", i, inOut[i], expected[i])
		}
	}
}

func TestClone(t *testing.T) {
	grid, _ := gridData.NewRGrid(-5.0, 5.0, 10)
	ke := NewKeDVR(grid, 2.0)

	ke.GetMat()
	err := ke.ensureEigen()
	if err != nil {
		return
	}

	keClone := ke.Clone()

	if keClone.Ndims() != ke.Ndims() {
		t.Errorf("Clone Ndims = %d, want %d", keClone.Ndims(), ke.Ndims())
	}

	if keClone.Mass() != ke.Mass() {
		t.Errorf("Clone Mass = %f, want %f", keClone.Mass(), ke.Mass())
	}

	ke.Clear()
	if keClone.eigVals == nil {
		t.Error("Clone should have independent eigen cache")
	}
}

func TestCloneWithoutCache(t *testing.T) {
	grid, _ := gridData.NewRGrid(-5.0, 5.0, 10)
	ke := NewKeDVR(grid, 2.0)
	keClone := ke.Clone()

	_, err, _ := keClone.RealDiagonalize()
	if err != nil {
		t.Fatalf("Clone RealDiagonalize failed: %v", err)
	}
}

func TestEigenvaluesOrdering(t *testing.T) {
	grid, _ := gridData.NewRGrid(-5.0, 5.0, 15)
	ke := NewKeDVR(grid, 1.0)

	eigenvalues, _, err := ke.RealDiagonalize()
	if err != nil {
		t.Fatalf("RealDiagonalize failed: %v", err)
	}

	// Check eigenvalues are in ascending order (gonum returns them sorted)
	for i := 1; i < len(eigenvalues); i++ {
		if eigenvalues[i] < eigenvalues[i-1]-tolerance {
			t.Errorf("Eigenvalues not sorted: eigenvalues[%d]=%f < eigenvalues[%d]=%f",
				i, eigenvalues[i], i-1, eigenvalues[i-1])
		}
	}
}

func TestExpDtConsistency(t *testing.T) {
	grid, _ := gridData.NewRGrid(-5.0, 5.0, 8)
	ke := NewKeDVR(grid, 1.0)

	in := make([]float64, 8)
	for i := range in {
		in[i] = math.Sin(float64(i) * 0.5)
	}

	// Apply exp(dt*K) twice with dt=0.1
	out1, _ := ke.ExpDt(0.1, in)
	out2, _ := ke.ExpDt(0.1, out1)

	// Apply exp(dt*K) once with dt=0.2
	out3, _ := ke.ExpDt(0.2, in)

	// Results should match: exp(0.2*K) = exp(0.1*K) * exp(0.1*K)
	for i := range out2 {
		if !floatEqual(out2[i], out3[i], 1e-8) {
			t.Errorf("Consistency check failed at [%d]: exp(0.1)^2 = %f, exp(0.2) = %f",
				i, out2[i], out3[i])
		}
	}
}

func TestLargeDimension(t *testing.T) {
	// Test with a larger grid to ensure no issues with larger matrices
	grid, _ := gridData.NewRGrid(-10.0, 10.0, 100)
	ke := NewKeDVR(grid, 1.0)

	if ke.Ndims() != 100 {
		t.Errorf("Ndims = %d, want 100", ke.Ndims())
	}

	// Should be able to compute matrix and eigenvalues
	eigenvalues, _, err := ke.RealDiagonalize()
	if err != nil {
		t.Fatalf("RealDiagonalize failed for large matrix: %v", err)
	}

	if len(eigenvalues) != 100 {
		t.Errorf("Expected 100 eigenvalues, got %d", len(eigenvalues))
	}
}

func TestDifferentMasses(t *testing.T) {
	grid, _ := gridData.NewRGrid(-5.0, 5.0, 10)

	// Heavier mass should give smaller kinetic energy eigenvalues
	keLightMass := NewKeDVR(grid, 1.0)
	keHeavyMass := NewKeDVR(grid, 10.0)

	evLight, _, _ := keLightMass.RealDiagonalize()
	evHeavy, _, _ := keHeavyMass.RealDiagonalize()

	// Each eigenvalue for heavy mass should be smaller (by factor of a mass ratio)
	for i := range evLight {
		ratio := evLight[i] / evHeavy[i]
		expectedRatio := 10.0
		if !floatEqual(ratio, expectedRatio, 0.01) {
			t.Errorf("Mass scaling incorrect at eigenvalue %d: ratio=%f, expected=%f",
				i, ratio, expectedRatio)
		}
	}
}

func BenchmarkGetMat(b *testing.B) {
	grid, _ := gridData.NewRGrid(-10.0, 10.0, 50)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		ke := NewKeDVR(grid, 1.0)
		ke.GetMat()
	}
}

func BenchmarkRealDiagonalize(b *testing.B) {
	grid, _ := gridData.NewRGrid(-10.0, 10.0, 50)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		ke := NewKeDVR(grid, 1.0)
		_, _, err := ke.RealDiagonalize()
		if err != nil {
			return
		}
	}
}

func BenchmarkExpDt(b *testing.B) {
	grid, _ := gridData.NewRGrid(-10.0, 10.0, 50)
	ke := NewKeDVR(grid, 1.0)

	in := make([]float64, 50)
	for i := range in {
		in[i] = math.Sin(float64(i) * 0.1)
	}

	err := ke.ensureEigen()
	if err != nil {
		return
	}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_, err := ke.ExpDt(0.01, in)
		if err != nil {
			return
		}
	}
}

func BenchmarkExpIdt(b *testing.B) {
	grid, _ := gridData.NewRGrid(-10.0, 10.0, 50)
	ke := NewKeDVR(grid, 1.0)

	in := make([]complex128, 50)
	for i := range in {
		in[i] = complex(math.Sin(float64(i)*0.1), math.Cos(float64(i)*0.1))
	}

	err := ke.ensureEigen()
	if err != nil {
		return
	}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_, err2 := ke.ExpIdt(0.01, in)
		if err2 != nil {
			return
		}
	}
}
