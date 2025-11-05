package NumericalRecipies

type RootFinder interface {
	Bisecion(a, b, tol float64) (float64, error)
	NewtonRaphson(x0, tol float64) (float64, error)
	secant(x0, x1, tol float64) (float64, error)
}