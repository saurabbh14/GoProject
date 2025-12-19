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
    // Single-interval / point methods
    Bisect(a, b, tol float64) (float64, error)
    Newton(x0, tol float64) (float64, error)
    Secant(x0, x1, tol float64) (float64, error)

    // Grid-level methods
    BisectionOnGrid(g *gridData.RadGrid, tol float64) ([]float64, error)
    NewtonRefineOnGrid(g *gridData.RadGrid, tol float64) ([]float64, error)
}

// rootFinderImpl implements RootFinder for a provided Rfunc.
type rootFinderImpl struct {
    f       gridData.Rfunc
    maxIter int
}

// NewRootFinder constructs a RootFinder.
func NewRootFinder(f gridData.Rfunc, maxIter int) RootFinder {
    if maxIter <= 0 {
        maxIter = 200
    }
    return &rootFinderImpl{f: f, maxIter: maxIter}
}

func (r *rootFinderImpl) Bisect(a, b, tol float64) (float64, error) {
    if a >= b {
        return 0, fmt.Errorf("invalid interval: a must be < b")
    }
    fa := r.f.EvaluateAt(a)
    fb := r.f.EvaluateAt(b)

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
    flo, fhi := fa, fb

    for i := 0; i < r.maxIter; i++ {
        m := 0.5 * (lo + hi)
        fm := r.f.EvaluateAt(m)

        if math.Abs(fm) <= tol || (hi-lo)/2 <= tol {
            return m, nil
        }

        if flo*fm <= 0 {
            hi = m
            fhi = fm
        } else {
            lo = m
            flo = fm
        }
    }

    return 0, fmt.Errorf("bisection did not converge after %d iterations", r.maxIter)
}

func (r *rootFinderImpl) Newton(x0, tol float64) (float64, error) {
    x := x0
    h := 1e-8
    for i := 0; i < r.maxIter; i++ {
        fx := r.f.EvaluateAt(x)
        if math.Abs(fx) <= tol {
            return x, nil
        }
        deriv := (r.f.EvaluateAt(x+h) - r.f.EvaluateAt(x-h)) / (2 * h)
        if deriv == 0 {
            return 0, errors.New("zero derivative encountered")
        }
        dx := fx / deriv
        x -= dx
        if math.Abs(dx) <= tol {
            return x, nil
        }
    }
    return 0, fmt.Errorf("Newton-Raphson did not converge after %d iterations", r.maxIter)
}

func (r *rootFinderImpl) Secant(x0, x1, tol float64) (float64, error) {
    xPrev := x0
    xCurr := x1
    fPrev := r.f.EvaluateAt(xPrev)
    fCurr := r.f.EvaluateAt(xCurr)

    for i := 0; i < r.maxIter; i++ {
        if math.Abs(fCurr) <= tol {
            return xCurr, nil
        }
        den := (fCurr - fPrev)
        if den == 0 {
            return 0, errors.New("division by zero in secant method")
        }
        xNext := xCurr - fCurr*(xCurr-xPrev)/den
        if math.Abs(xNext-xCurr) <= tol {
            return xNext, nil
        }
        xPrev, fPrev = xCurr, fCurr
        xCurr = xNext
        fCurr = r.f.EvaluateAt(xCurr)
    }
    return 0, fmt.Errorf("secant did not converge after %d iterations", r.maxIter)
}

func (r *rootFinderImpl) BisectionOnGrid(g *gridData.RadGrid, tol float64) ([]float64, error) {
    xs := g.RValues()
    if len(xs) < 2 {
        return nil, errors.New("grid must contain at least two points")
    }

    var roots []float64
    for i := 0; i < len(xs)-1; i++ {
        a := xs[i]
        b := xs[i+1]
        fa := r.f.EvaluateAt(a)
        fb := r.f.EvaluateAt(b)

        // exact root at grid point
        if math.Abs(fa) == 0 {
            roots = append(roots, a)
            continue
        }

        // sign change indicates a root in (a,b)
        if fa*fb < 0 {
            root, err := r.Bisect(a, b, tol)
            if err == nil {
                roots = append(roots, root)
            }
            // if bisect fails for a bracket, skip it
        }
    }

    return roots, nil
}

func (r *rootFinderImpl) NewtonRefineOnGrid(g *gridData.RadGrid, tol float64) ([]float64, error) {
    roots, err := r.BisectionOnGrid(g, tol)
    if err != nil {
        return nil, err
    }

    var refined []float64
    for _, x0 := range roots {
        rr, err := r.Newton(x0, tol)
        if err == nil {
            refined = append(refined, rr)
        } else {
            // if Newton fails, keep bisection result
            refined = append(refined, x0)
        }
    }
    return refined, nil
}