package NumericalRecipies

import (
    "errors"
    "fmt"
    "math"

    "GoProject/gridData"
)

// RootFinder is an interface that can find roots for a single interval
// (single-point routines) and over an entire RadGrid (grid-level routines).
type RootFinder interface {
    FindRoot(a, b, tol float64, g gridData.Rfunc) (float64, error)
    FindRootInPlace()
    FindRootsOnGrid()
}

// Bisection Struct for root finding
type Bisection struct {
    maxIter int 
}

func NewBisection(maxIter int) *Bisection {
    if maxIter <= 0 {
        maxIter = 200
    }
    return &Bisection{maxIter: maxIter}
}

func (r *Bisection) FindRoot(a, b, tol float64, f gridData.Rfunc) (float64, error) {
    if a >= b {
        return 0, fmt.Errorf("invalid interval: a must be < b")
    }
    fa := f.EvaluateAt(a)
    fb := f.EvaluateAt(b)

    if math.Abs(fa) == 0 {
        return a, nil
    }
    if math.Abs(fb) == 0 {
        return b, nil
    }
    if fa*fb > 0 {
        return 0, errors.New("endpoints do not bracket a root")
    }

    lo, hi := a, b
    flo := fa

    for i := 0; i < r.maxIter; i++ {
        m := 0.5 * (lo + hi)
        fm := f.EvaluateAt(m)

        if math.Abs(fm) <= tol || (hi-lo)/2 <= tol {
            return m, nil
        }

        if flo*fm <= 0 {
            hi = m
        } else {
            lo = m
            flo = fm
        }
    }

    return 0, fmt.Errorf("bisection did not converge after %d iterations", r.maxIter)
}