package OperatorAlgebra

import (
	"GoProject/gridData"
	"math"
	"testing"

	"gonum.org/v1/gonum/mat"
)

// TestNewFiniteDiff tests the constructor
func TestNewFiniteDiff(t *testing.T) {
	tests := []struct {
		name      string
		mass      float64
		order     int
		wantError bool
	}{
		{"Valid order 3", 1.0, 3, false},
		{"Valid order 5", 1.0, 5, false},
		{"Valid order 6", 1.0, 6, false},
		{"Invalid order 2", 1.0, 2, true},
		{"Invalid order 4", 1.0, 4, true},
		{"Invalid order 7", 1.0, 7, true},
	}

	grid, err := gridData.NewFromLength(20, 100)
	if err != nil {
		t.Errorf("NewFromLength failed with error: %v", err)
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			fd, err := NewFiniteDiff(grid, tt.mass, tt.order)

			if tt.wantError {
				if err == nil {
					t.Errorf("NewFiniteDiff() expected error but got none")
				}
				if fd != nil {
					t.Errorf("NewFiniteDiff() expected nil FiniteDiff on error")
				}
			} else {
				if err != nil {
					t.Errorf("NewFiniteDiff() unexpected error: %v", err)
				}
				if fd == nil {
					t.Errorf("NewFiniteDiff() expected non-nil FiniteDiff")
				}
			}
		})
	}
}

// TestKineticCoefficient verifies the kinetic coefficient calculation
func TestKineticCoefficient(t *testing.T) {
	grid, err := gridData.NewFromLength(20, 100)
	if err != nil {
		t.Errorf("NewFromLength failed with error: %v", err)
	}
	mass := 2.0

	fd, err := NewFiniteDiff(grid, mass, 3)
	if err != nil {
		t.Fatalf("NewFiniteDiff() error: %v", err)
	}

	dxSq := grid.DeltaR() * grid.DeltaR()
	expected := -1.0 / (2.0 * mass * dxSq)

	if math.Abs(fd.kineticCoeff-expected) > 1e-10 {
		t.Errorf("kineticCoeff = %v, want %v", fd.kineticCoeff, expected)
	}
}

// TestCompute3rdOrderCoefficients tests 3-point stencil
func TestCompute3rdOrderCoefficients(t *testing.T) {
	grid, err := gridData.NewFromLength(20, 100)
	if err != nil {
		t.Errorf("NewFromLength failed with error: %v", err)
	}
	fd, _ := NewFiniteDiff(grid, 1.0, 3)

	keMat := mat.NewDense(5, 5, nil)
	fd.compute3rdOrderCoefficients(keMat)

	expectedDiag := -2.0 * fd.kineticCoeff
	for i := 0; i < 5; i++ {
		if math.Abs(keMat.At(i, i)-expectedDiag) > 1e-10 {
			t.Errorf("Diagonal[%d] = %v, want %v", i, keMat.At(i, i), expectedDiag)
		}
	}

	for i := 0; i < 4; i++ {
		if math.Abs(keMat.At(i, i+1)-fd.kineticCoeff) > 1e-10 {
			t.Errorf("Upper diagonal[%d] = %v, want %v", i, keMat.At(i, i+1), fd.kineticCoeff)
		}
		if math.Abs(keMat.At(i+1, i)-fd.kineticCoeff) > 1e-10 {
			t.Errorf("Lower diagonal[%d] = %v, want %v", i, keMat.At(i+1, i), fd.kineticCoeff)
		}
	}

	if !isSymmetric(keMat) {
		t.Error("Matrix is not symmetric")
	}
}

// TestCompute5thOrderCoefficients tests 5-point stencil
func TestCompute5thOrderCoefficients(t *testing.T) {
	grid, err := gridData.NewFromLength(20, 100)
	if err != nil {
		t.Errorf("NewFromLength failed with error: %v", err)
	}
	fd, _ := NewFiniteDiff(grid, 1.0, 5)

	keMat := mat.NewDense(6, 6, nil)
	fd.compute5thOrderCoefficients(keMat)

	k := fd.kineticCoeff

	for i := 0; i < 6; i++ {
		expected := 2.5 * k
		if math.Abs(keMat.At(i, i)-expected) > 1e-10 {
			t.Errorf("Diagonal[%d] = %v, want %v", i, keMat.At(i, i), expected)
		}
	}

	for i := 0; i < 5; i++ {
		expected := -4.0 / 3.0 * k
		if math.Abs(keMat.At(i, i+1)-expected) > 1e-10 {
			t.Errorf("First off-diag[%d] = %v, want %v", i, keMat.At(i, i+1), expected)
		}
	}

	for i := 0; i < 4; i++ {
		expected := 1.0 / 12.0 * k
		if math.Abs(keMat.At(i, i+2)-expected) > 1e-10 {
			t.Errorf("Second off-diag[%d] = %v, want %v", i, keMat.At(i, i+2), expected)
		}
	}

	if !isSymmetric(keMat) {
		t.Error("Matrix is not symmetric")
	}
}

// TestCompute7thOrderCoefficients tests 7-point stencil
func TestCompute7thOrderCoefficients(t *testing.T) {
	grid, err := gridData.NewFromLength(20, 100)
	if err != nil {
		t.Errorf("NewFromLength failed with error: %v", err)
	}
	fd, _ := NewFiniteDiff(grid, 1.0, 6) // Note: using order 6 to call 7th order

	keMat := mat.NewDense(8, 8, nil)
	fd.compute7thOrderCoefficients(keMat)

	k := fd.kineticCoeff
	expectedDiag := 49.0 / 18.0 * k
	if math.Abs(keMat.At(0, 0)-expectedDiag) > 1e-10 {
		t.Errorf("Diagonal = %v, want %v", keMat.At(0, 0), expectedDiag)
	}

	expectedOff1 := -1.5 * k
	if math.Abs(keMat.At(0, 1)-expectedOff1) > 1e-10 {
		t.Errorf("First off-diag = %v, want %v", keMat.At(0, 1), expectedOff1)
	}

	expectedOff2 := 3.0 / 20.0 * k
	if math.Abs(keMat.At(0, 2)-expectedOff2) > 1e-10 {
		t.Errorf("Second off-diag = %v, want %v", keMat.At(0, 2), expectedOff2)
	}

	expectedOff3 := -1.0 / 90.0 * k
	if math.Abs(keMat.At(0, 3)-expectedOff3) > 1e-10 {
		t.Errorf("Third off-diag = %v, want %v", keMat.At(0, 3), expectedOff3)
	}

	if !isSymmetric(keMat) {
		t.Error("Matrix is not symmetric")
	}
}

// TestCompute9thOrderCoefficients tests 9-point stencil
func TestCompute9thOrderCoefficients(t *testing.T) {
	grid, err := gridData.NewFromLength(20, 100)
	if err != nil {
		t.Errorf("NewFromLength failed with error: %v", err)
	}
	fd, _ := NewFiniteDiff(grid, 1.0, 5)

	keMat := mat.NewDense(10, 10, nil)
	fd.compute9thOrderCoefficients(keMat)

	k := fd.kineticCoeff
	coeffs := []struct {
		offset int
		value  float64
	}{
		{0, 205.0 / 72.0},
		{1, -8.0 / 5.0},
		{2, 1.0 / 5.0},
		{3, -8.0 / 315.0},
		{4, 1.0 / 560.0},
	}

	for _, c := range coeffs {
		expected := c.value * k
		for i := 0; i < 10-c.offset; i++ {
			if math.Abs(keMat.At(i, i+c.offset)-expected) > 1e-10 {
				t.Errorf("Offset %d at [%d,%d] = %v, want %v",
					c.offset, i, i+c.offset, keMat.At(i, i+c.offset), expected)
			}
		}
	}

	if !isSymmetric(keMat) {
		t.Error("Matrix is not symmetric")
	}
}

// TestComputeMatrix tests the dispatcher method
func TestComputeMatrix(t *testing.T) {
	grid, err := gridData.NewFromLength(20, 100)
	if err != nil {
		t.Errorf("NewFromLength failed with error: %v", err)
	}

	orders := []int{3, 5}
	for _, order := range orders {
		t.Run("Order"+string(rune(order+'0')), func(t *testing.T) {
			fd, err := NewFiniteDiff(grid, 1.0, order)
			if err != nil {
				t.Fatalf("NewFiniteDiff() error: %v", err)
			}

			keMat := mat.NewDense(5, 5, nil)
			fd.computeMatrix(keMat)

			hasNonZero := false
			for i := 0; i < 5; i++ {
				for j := 0; j < 5; j++ {
					if math.Abs(keMat.At(i, j)) > 1e-10 {
						hasNonZero = true
						break
					}
				}
			}

			if !hasNonZero {
				t.Error("Matrix is all zeros")
			}
		})
	}
}

// Helper function to check matrix symmetry
func isSymmetric(m *mat.Dense) bool {
	r, c := m.Dims()
	if r != c {
		return false
	}

	for i := 0; i < r; i++ {
		for j := i + 1; j < c; j++ {
			if math.Abs(m.At(i, j)-m.At(j, i)) > 1e-10 {
				return false
			}
		}
	}
	return true
}

// Benchmark tests
func BenchmarkCompute3rdOrder(b *testing.B) {
	grid, err := gridData.NewFromLength(20, 100)
	if err != nil {
		b.Errorf("NewFromLength failed with error: %v", err)
	}

	fd, _ := NewFiniteDiff(grid, 1.0, 3)
	keMat := mat.NewDense(100, 100, nil)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		fd.compute3rdOrderCoefficients(keMat)
	}
}

func BenchmarkCompute9thOrder(b *testing.B) {
	grid, err := gridData.NewFromLength(20, 100)
	if err != nil {
		b.Errorf("NewFromLength failed with error: %v", err)
	}

	fd, _ := NewFiniteDiff(grid, 1.0, 5)
	keMat := mat.NewDense(100, 100, nil)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		fd.compute9thOrderCoefficients(keMat)
	}
}
