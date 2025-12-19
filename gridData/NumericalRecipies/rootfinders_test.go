package NumericalRecipies

import (
    "math"
    "testing"

    "GoProject/gridData"
)

func approxEqual(a, b, tol float64) bool { return math.Abs(a-b) <= tol }

func TestBisect_Newton_Secant_Sqrt2(t *testing.T) {
    f := gridData.RFuncFunc(func(x float64) float64 { return x*x - 2 })
    rf := NewRootFinder(f, 200)

    // Bisect
    r, err := rf.Bisect(1, 2, 1e-12)
    if err != nil {
        t.Fatalf("Bisect error: %v", err)
    }
    if !approxEqual(r, math.Sqrt(2), 1e-9) {
        t.Fatalf("Bisect root mismatch: got %v want %v", r, math.Sqrt(2))
    }

    // Newton
    rn, err := rf.Newton(1.5, 1e-12)
    if err != nil {
        t.Fatalf("Newton error: %v", err)
    }
    if !approxEqual(rn, math.Sqrt(2), 1e-9) {
        t.Fatalf("Newton root mismatch: got %v want %v", rn, math.Sqrt(2))
    }

    // Secant
    rs, err := rf.Secant(1.0, 2.0, 1e-12)
    if err != nil {
        t.Fatalf("Secant error: %v", err)
    }
    if !approxEqual(rs, math.Sqrt(2), 1e-9) {
        t.Fatalf("Secant root mismatch: got %v want %v", rs, math.Sqrt(2))
    }
}

func TestBisect_EndpointAndNoBracket(t *testing.T) {
    // endpoint root should be returned when an endpoint is exact root
    f := gridData.RFuncFunc(func(x float64) float64 { return x })
    rf := NewRootFinder(f, 100)
    r, err := rf.Bisect(-1, 0, 1e-12)
    if err != nil {
        t.Fatalf("expected nil err for endpoint root, got %v", err)
    }
    if !approxEqual(r, 0, 1e-12) {
        t.Fatalf("expected root 0 got %v", r)
    }

    // non-bracketing interval must return an error
    g := gridData.RFuncFunc(func(x float64) float64 { return x*x + 1 })
    rf2 := NewRootFinder(g, 100)
    _, err = rf2.Bisect(-1, 1, 1e-12)
    if err == nil {
        t.Fatalf("expected error for non-bracketing interval")
    }
}

func TestBisectionOnGrid_SinRoots(t *testing.T) {
    g, err := gridData.NewRGrid(-4*math.Pi, 4*math.Pi, 2000)
    if err != nil {
        t.Fatalf("failed to create grid: %v", err)
    }

    f := gridData.RFuncFunc(func(x float64) float64 { return math.Sin(x) })
    rf := NewRootFinder(f, 200)
    roots, err := rf.BisectionOnGrid(g, 1e-9)
    if err != nil {
        t.Fatalf("BisectionOnGrid error: %v", err)
    }

    // expect 9 roots at k*pi for k=-4..4
    if len(roots) < 9 {
        t.Fatalf("expected at least 9 roots, got %d", len(roots))
    }

    for _, r := range roots {
        if math.Abs(math.Sin(r)) > 1e-6 {
            t.Fatalf("root not accurate: sin(%v) = %g", r, math.Sin(r))
        }
    }
}

func TestNewtonRefineOnGrid(t *testing.T) {
    g, err := gridData.NewRGrid(-4*math.Pi, 4*math.Pi, 2000)
    if err != nil {
        t.Fatalf("failed to create grid: %v", err)
    }

    f := gridData.RFuncFunc(func(x float64) float64 { return math.Sin(x) })
    rf := NewRootFinder(f, 200)
    refined, err := rf.NewtonRefineOnGrid(g, 1e-12)
    if err != nil {
        t.Fatalf("NewtonRefineOnGrid error: %v", err)
    }

    if len(refined) < 9 {
        t.Fatalf("expected at least 9 refined roots, got %d", len(refined))
    }

    for _, r := range refined {
        if math.Abs(math.Sin(r)) > 1e-9 {
            t.Fatalf("refined root not accurate enough: sin(%v) = %g", r, math.Sin(r))
        }
    }
}

func TestSecant_SameInitialsError(t *testing.T) {
    f := gridData.RFuncFunc(func(x float64) float64 { return x*x - 2 })
    rf := NewRootFinder(f, 100)
    _, err := rf.Secant(1, 1, 1e-12)
    if err == nil {
        t.Fatalf("expected error for secant with identical initial guesses")
    }
}