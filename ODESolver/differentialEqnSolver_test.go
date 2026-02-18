package EquationSolver

import (
	"fmt"
	"math"
	"testing"
)

// a simple decay equation: dx/dt = -x
type MockPotential struct{}

// EvaluateOnRGridTime implements [gridData.TDPotentialOp].
func (m *MockPotential) EvaluateOnRGridTime(x []float64, t float64) []float64 {
	panic("unimplemented")
}

func (m *MockPotential) EvaluateAtTime(x, t float64) float64 {
	// dx/dt = -x
	return -x
}

func (m *MockPotential) EvaluateOnRGridTimeInPlace(x []float64, buff []float64, t float64) {
	for i, v := range x {
		buff[i] = -v
	}
}

// TestEulerExplicit_NextStep tests the scalar stepping function.
func TestEulerExplicit_NextStep(t *testing.T) {
	mockFunc := &MockPotential{}

	// Problem: dx/dt = -x, x(0) = 1.0
	// Euler explicit: x_new = x + dt * (-x) = x * (1 - dt)
	initialX := 1.0
	dt := 0.1
	expectedX := initialX * (1.0 - dt)

	solver := &EulerExplicit{}
	solver = solver.NewDefine(dt, mockFunc)

	// Step 1
	nextX, err := solver.NextStep(initialX, 0.0)
	if err != nil {
		t.Fatalf("NextStep failed: %v", err)
	}

	if math.Abs(nextX-expectedX) > 1e-9 {
		t.Errorf("NextStep wrong result. Got %f, expected %f", nextX, expectedX)
	}

	// Step 2
	expectedX2 := expectedX * (1.0 - dt)
	nextX2, err := solver.NextStep(nextX, dt)
	if err != nil {
		t.Fatalf("NextStep failed on second step: %v", err)
	}

	if math.Abs(nextX2-expectedX2) > 1e-9 {
		t.Errorf("NextStep wrong result on step 2. Got %f, expected %f", nextX2, expectedX2)
	}
}

// TestEulerExplicit_NextStepOnGrid tests the vector/grid stepping function.
func TestEulerExplicit_NextStepOnGrid(t *testing.T) {
	mockFunc := &MockPotential{}
	dt := 0.1
	solver := &EulerExplicit{}
	solver = solver.NewDefine(dt, mockFunc)

	// Grid inputs
	gridX := []float64{1.0, 2.0, 3.0}
	// Expected: x_new = x * (1 - dt)
	expectedGrid := []float64{
		1.0 * (1.0 - dt),
		2.0 * (1.0 - dt),
		3.0 * (1.0 - dt),
	}

	err := solver.NextStepOnGrid(gridX, 0.0)
	if err != nil {
		t.Fatalf("NextStepOnGrid failed: %v", err)
	}
	fmt.Println(gridX)

	for i, val := range gridX {
		if math.Abs(val-expectedGrid[i]) > 1e-9 {
			t.Errorf("Index %d: got %f, expected %f", i, val, expectedGrid[i])
		}
	}
}

// TestEulerExplicit_100StepsStable checks long-run behavior in a stable regime (dt < 2).
func TestEulerExplicit_100StepsStable(t *testing.T) {
    mockFunc := &MockPotential{}
    dt := 0.01
    steps := 100

    solver := (&EulerExplicit{}).NewDefine(dt, mockFunc)

    x := 1.0
    time := 0.0
    for i := 0; i < steps; i++ {
        var err error
        x, err = solver.NextStep(x, time)
        if err != nil {
            t.Fatalf("NextStep failed at step %d: %v", i, err)
        }
        time += dt
    }

    // exact Euler recurrence for this ODE: x_n = (1-dt)^n x0
    expectedEuler := math.Pow(1.0-dt, float64(steps))
    if math.Abs(x-expectedEuler) > 1e-12 {
        t.Errorf("100-step mismatch against Euler closed form. got=%g expected=%g", x, expectedEuler)
    }

    // optional sanity against analytical x(t)=exp(-t), allow looser tolerance (method error)
    expectedAnalytical := math.Exp(-float64(steps) * dt)
    if math.Abs(x-expectedAnalytical) > 2e-3 {
        t.Errorf("100-step mismatch against analytical solution. got=%g expected=%g", x, expectedAnalytical)
    }
}

// TestEulerExplicit_100StepsDivergence checks divergence in unstable regime (dt > 2).
func TestEulerExplicit_100StepsDivergence(t *testing.T) {
    mockFunc := &MockPotential{}
    dt := 2.1 // amplification factor = (1-dt) = -1.1, |factor| > 1 => divergence
    steps := 100

    solver := (&EulerExplicit{}).NewDefine(dt, mockFunc)

    x := 1.0
    time := 0.0
    for i := 0; i < steps; i++ {
        var err error
        x, err = solver.NextStep(x, time)
        if err != nil {
            t.Fatalf("NextStep failed at step %d: %v", i, err)
        }
        time += dt
    }

    expected := math.Pow(1.0-dt, float64(steps))
    if math.Abs((x-expected)/expected) > 1e-12 {
        t.Errorf("divergence recurrence mismatch. got=%g expected=%g", x, expected)
    }

    if math.Abs(x) <= 1e3 {
        t.Errorf("expected divergence after %d steps, got |x|=%g", steps, math.Abs(x))
    }
}

// TestEulerExplicit_ReDefine verifies that changing parameters works.
func TestEulerExplicit_ReDefine(t *testing.T) {
	mockFunc := &MockPotential{}
	dt1 := 0.1
	solver := &EulerExplicit{}
	solver = solver.NewDefine(dt1, mockFunc)

	// Check initial set up
	if solver.deltaTime != dt1 {
		t.Errorf("Initial dt mismatch. Got %f, expected %f", solver.deltaTime, dt1)
	}

	dt2 := 0.05
	solver.ReDefine(dt2, mockFunc)

	if solver.deltaTime != dt2 {
		t.Errorf("Redefined dt mismatch. Got %f, expected %f", solver.deltaTime, dt2)
	}
}

// TestEulerExplicit_Name verifies the Name method.
func TestEulerExplicit_Name(t *testing.T) {
	solver := &EulerExplicit{}
	expected := "Euler Explicit Method"
	if solver.Name() != expected {
		t.Errorf("Name mismatch. Got %s, expected %s", solver.Name(), expected)
	}
}