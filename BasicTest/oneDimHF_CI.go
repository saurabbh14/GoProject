package BasicTest

import (
	"GoProject/gridData"
	"fmt"
	"math"

	"gonum.org/v1/gonum/mat"
)

type SineBasisHF struct {
	grid *gridData.RadGrid
	poten gridData.PotentialOp[float64]
	ne   uint8
	nAos int
}

func NewSinBasis(grid *gridData.RadGrid, pot gridData.PotentialOp[float64], ne uint8, nAos int) *SineBasisHF {
	return &SineBasisHF{
		grid: grid,
		poten: pot,
		ne:   ne,
		nAos: nAos,
	}
}

// Basis phi_n(x) = sqrt(2/L) * sin(n*pi*(x+L/2)/L), n=1..M
func (s *SineBasisHF) Basis() *mat.Dense {
	Ngrid := int(s.grid.NPoints())
	Norm := math.Sqrt(2.0/s.grid.Length())

	Phi := mat.NewDense(s.nAos, Ngrid, nil)
	xp := s.grid.RValues()
	for iAO := 0; iAO < s.nAos; iAO++ {
		nPibyL := (float64(iAO+ 1)*math.Pi)/s.grid.Length()
		for i := 0; i < Ngrid; i++ {
			Phi.Set(iAO, i, Norm * math.Sin(nPibyL*xp[i]))
		}
	}
	return Phi
}

func (s *SineBasisHF) BuildHcoreIntegrals(){
	Phi := s.Basis()
	v := s.grid.PotentialOnGrid(s.poten)
	piByL := math.Pi/s.grid.Length()

	h := mat.NewDense(s.nAos, s.nAos, nil)
	for n := 0; n < s.nAos; n++ {
		nPibyL := float64(n+ 1)*piByL
		h.Set(n, n, 0.5 * math.Pow(nPibyL, 2))
	}

	V := mat.NewDense(s.nAos, s.nAos, nil)
	for p := 0; p < s.nAos; p++ {
		for q := 0; q < s.nAos; q++ {
			sum := 0.0
			for i := 0; i < int(s.grid.NPoints()); i++ {
				sum += Phi.At(p, i) * v[i] * Phi.At(q, i) * s.grid.DeltaR()
			}
			V.Set(p, q, sum)
		}
	}
	h.Add(h, V)
}

// BuildIntegrals constructs one-electron and two-electron integrals
func (s *SineBasisHF) BuildIntegrals() ([]float64, *mat.Dense, [][][][]float64) {
	nGrid := int(s.grid.NPoints())
	W := mat.NewDense(nGrid, nGrid, nil)
	for i := 0; i < nGrid; i++ {
		for j := 0; j < nGrid; j++ {
			diff := x[i] - x[j]
			W.Set(i, j, strength/math.Sqrt(diff*diff+aSoft*aSoft))
		}
	}

	// Build g[p,q,r,s] tensor
	// A[p,q,i] = phi_p(i) * phi_q(i)
	A := make([][][]float64, s.nAos)
	for p := 0; p < s.nAos; p++ {
		A[p] = make([][]float64, s.nAos)
		for q := 0; q < M; q++ {
			A[p][q] = make([]float64, Ngrid)
			for i := 0; i < Ngrid; i++ {
				A[p][q][i] = Phi.At(p, i) * Phi.At(q, i)
			}
		}
	}

	// B[p,q,j] = sum_i A[p,q,i] * W[i,j] * dx
	B := make([][][]float64, M)
	for p := 0; p < M; p++ {
		B[p] = make([][]float64, M)
		for q := 0; q < M; q++ {
			B[p][q] = make([]float64, Ngrid)
			for j := 0; j < Ngrid; j++ {
				sum := 0.0
				for i := 0; i < Ngrid; i++ {
					sum += A[p][q][i] * W.At(i, j) * dx
				}
				B[p][q][j] = sum
			}
		}
	}

	// g[p,q,r,s] = sum_j B[p,q,j] * A[r,s,j] * dx
	g := make([][][][]float64, M)
	for p := 0; p < M; p++ {
		g[p] = make([][][]float64, M)
		for q := 0; q < M; q++ {
			g[p][q] = make([][]float64, M)
			for r := 0; r < M; r++ {
				g[p][q][r] = make([]float64, M)
				for s := 0; s < M; s++ {
					sum := 0.0
					for j := 0; j < Ngrid; j++ {
						sum += B[p][q][j] * A[r][s][j] * dx
					}
					g[p][q][r][s] = sum
				}
			}
		}
	}

	return x, h, g
}

// RHF performs restricted Hartree-Fock calculation
func (s *SineBasisHF) RHF(h *mat.Dense, g [][][][]float64, Nelec int, maxIter int, conv, damp float64, verbose bool) (float64, []float64, *mat.Dense, *mat.Dense, *mat.Dense) {
	if Nelec%2 != 0 {
		panic("RHF assumes closed shell: Nelec must be even")
	}

	nocc := Nelec / 2
	M := h.RawMatrix().Rows

	// Convert Dense to Symmetric
	sym := mat.NewSymDense(h., nil)
	for i := 0; i < k.ndims; i++ {
		for j := i; j < k.ndims; j++ {
			sym.SetSym(i, j, KM.At(i, j))
		}
	}

	// Initial guess: diagonalize h
	var eig mat.EigenSym
	ok := eig.Factorize(h, true)
	if !ok {
		panic("eigenvalue decomposition failed")
	}

	eps := make([]float64, M)
	eig.Values(eps)

	C := mat.NewDense(M, M, nil)
	eig.VectorsTo(C)

	// Density matrix: P = 2 * Cocc * Cocc^T
	Cocc := mat.NewDense(M, nocc, nil)
	Cocc.Copy(C.Slice(0, M, 0, nocc))

	P := mat.NewDense(M, M, nil)
	P.Product(Cocc, Cocc.T())
	P.Scale(2.0, P)

	var EOld *float64

	for it := 1; it <= maxIter; it++ {
		// Coulomb: J_pq = sum_rs P_rs * g[p,q,r,s]
		J := mat.NewDense(M, M, nil)
		for p := 0; p < M; p++ {
			for q := 0; q < M; q++ {
				sum := 0.0
				for r := 0; r < M; r++ {
					for s := 0; s < M; s++ {
						sum += P.At(r, s) * g[p][q][r][s]
					}
				}
				J.Set(p, q, sum)
			}
		}

		// Exchange: K_pq = sum_rs P_rs * g[p,r,q,s]
		K := mat.NewDense(M, M, nil)
		for p := 0; p < M; p++ {
			for q := 0; q < M; q++ {
				sum := 0.0
				for r := 0; r < M; r++ {
					for s := 0; s < M; s++ {
						sum += P.At(r, s) * g[p][r][q][s]
					}
				}
				K.Set(p, q, sum)
			}
		}

		// Fock matrix: F = h + J - 0.5*K
		F := mat.NewDense(M, M, nil)
		F.Copy(h)
		F.Add(F, J)
		K.Scale(0.5, K)
		F.Sub(F, K)

		// Diagonalize Fock matrix
		var eigF mat.EigenSym
		ok := eigF.Factorize(F, true)
		if !ok {
			panic("Fock matrix diagonalization failed")
		}

		eigF.Values(eps)

		CNew := mat.NewDense(M, M, nil)
		eigF.VectorsTo(CNew)

		// New density matrix
		CoccNew := mat.NewDense(M, nocc, nil)
		CoccNew.Copy(CNew.Slice(0, M, 0, nocc))

		PNew := mat.NewDense(M, M, nil)
		PNew.Product(CoccNew, CoccNew.T())
		PNew.Scale(2.0, PNew)

		// Damping
		PDamped := mat.NewDense(M, M, nil)
		PDamped.Scale(1.0-damp, PNew)
		PTemp := mat.NewDense(M, M, nil)
		PTemp.Scale(damp, P)
		PDamped.Add(PDamped, PTemp)
		P.Copy(PDamped)

		// Energy: E = sum(P * h) + 0.5 * sum(P * G), where G = F - h
		G := mat.NewDense(M, M, nil)
		G.Sub(F, h)

		E := 0.0
		for i := 0; i < M; i++ {
			for j := 0; j < M; j++ {
				E += P.At(i, j) * h.At(i, j)
				E += 0.5 * P.At(i, j) * G.At(i, j)
			}
		}

		// Convergence check
		dP := mat.NewDense(M, M, nil)
		dP.Sub(PNew, P)
		dPNorm := mat.Norm(dP, 2)

		var dE *float64
		if EOld != nil {
			diff := math.Abs(E - *EOld)
			dE = &diff
		}

		if verbose {
			if dE == nil {
				fmt.Printf("RHF iter %3d: E = %.12f\n", it, E)
			} else {
				fmt.Printf("RHF iter %3d: E = %.12f  |dE|=%.3e  ||dP||=%.3e\n", it, E, *dE, dPNorm)
			}
		}

		if dE != nil && *dE < conv && dPNorm < 1e-8 {
			return E, eps, CNew, PNew, F
		}

		EOld = &E
		C.Copy(CNew)
	}

	panic("RHF did not converge. Try larger damp or better guess.")
}

// ----------------------------
// FCI (Full Configuration Interaction)
// ----------------------------

func popcount(x uint64) int {
	count := 0
	for x != 0 {
		count++
		x &= x - 1
	}
	return count
}

func occIndices(det uint64, nspin int) []int {
	occ := make([]int, 0)
	for i := 0; i < nspin; i++ {
		if (det>>uint(i))&1 == 1 {
			occ = append(occ, i)
		}
	}
	return occ
}

func phaseSingle(det uint64, i, a int) int {
	if i == a {
		return 1
	}

	var mask uint64
	if i < a {
		mask = ((1 << uint(a)) - 1) ^ ((1 << uint(i+1)) - 1)
	} else {
		mask = ((1 << uint(i)) - 1) ^ ((1 << uint(a+1)) - 1)
	}

	n := popcount(det & mask)
	if n%2 == 0 {
		return 1
	}
	return -1
}

func phaseDouble(det uint64, i, j, a, b int) int {
	sign := 1
	d := det

	// First excitation i->a
	sign *= phaseSingle(d, i, a)
	d = d ^ (1 << uint(i)) ^ (1 << uint(a))

	// Second excitation j->b
	sign *= phaseSingle(d, j, b)
	return sign
}

func spatialToSpinOrbIntegrals(h *mat.Dense, g [][][][]float64) (*mat.Dense, [][][][]float64) {
	M := h.RawMatrix().Rows
	nspin := 2 * M

	// Spin-orbital one-electron integrals
	hso := mat.NewDense(nspin, nspin, nil)
	for p := 0; p < nspin; p++ {
		var spP int
		if p < M {
			spP = 0
		} else {
			spP = 1
		}
		mu := p
		if p >= M {
			mu = p - M
		}

		for q := 0; q < nspin; q++ {
			var spQ int
			if q < M {
				spQ = 0
			} else {
				spQ = 1
			}
			nu := q
			if q >= M {
				nu = q - M
			}

			if spP == spQ {
				hso.Set(p, q, h.At(mu, nu))
			}
		}
	}

	// Antisymmetrized two-electron integrals
	eri := make([][][][]float64, nspin)
	for p := 0; p < nspin; p++ {
		eri[p] = make([][][]float64, nspin)
		for q := 0; q < nspin; q++ {
			eri[p][q] = make([][]float64, nspin)
			for r := 0; r < nspin; r++ {
				eri[p][q][r] = make([]float64, nspin)
			}
		}
	}

	for p := 0; p < nspin; p++ {
		var spP int
		if p < M {
			spP = 0
		} else {
			spP = 1
		}
		mu := p
		if p >= M {
			mu = p - M
		}

		for q := 0; q < nspin; q++ {
			var spQ int
			if q < M {
				spQ = 0
			} else {
				spQ = 1
			}
			nu := q
			if q >= M {
				nu = q - M
			}

			for r := 0; r < nspin; r++ {
				var spR int
				if r < M {
					spR = 0
				} else {
					spR = 1
				}
				lam := r
				if r >= M {
					lam = r - M
				}

				for s := 0; s < nspin; s++ {
					var spS int
					if s < M {
						spS = 0
					} else {
						spS = 1
					}
					sig := s
					if s >= M {
						sig = s - M
					}

					v1 := 0.0
					if spP == spR && spQ == spS {
						v1 = g[mu][lam][nu][sig]
					}

					v2 := 0.0
					if spP == spS && spQ == spR {
						v2 = g[mu][sig][nu][lam]
					}

					eri[p][q][r][s] = v1 - v2
				}
			}
		}
	}

	return hso, eri
}

func combinations(n, k int) [][]int {
	result := [][]int{}

	var generate func(start int, combo []int)
	generate = func(start int, combo []int) {
		if len(combo) == k {
			temp := make([]int, k)
			copy(temp, combo)
			result = append(result, temp)
			return
		}

		for i := start; i < n; i++ {
			generate(i+1, append(combo, i))
		}
	}

	generate(0, []int{})
	return result
}

func determinants(nspin, Nelec int) []uint64 {
	dets := []uint64{}
	combs := combinations(nspin, Nelec)

	for _, occ := range combs {
		var d uint64
		for _, i := range occ {
			d |= 1 << uint(i)
		}
		dets = append(dets, d)
	}

	return dets
}

func FCI(h *mat.Dense, g [][][][]float64, Nelec int, verbose bool) (float64, []float64, *mat.Dense, []uint64) {
	hso, eri := spatialToSpinOrbIntegrals(h, g)
	nspin := hso.RawMatrix().Rows

	dets := determinants(nspin, Nelec)
	ndet := len(dets)

	if verbose {
		fmt.Printf("FCI space: nspin=%d, Nelec=%d, ndet=%d\n", nspin, Nelec, ndet)
	}

	// Create index map
	idx := make(map[uint64]int)
	for i, d := range dets {
		idx[d] = i
	}

	// Build Hamiltonian matrix
	H := mat.NewDense(ndet, ndet, nil)

	// Diagonal elements
	for I, d := range dets {
		occ := occIndices(d, nspin)

		// One-electron diagonal
		e1 := 0.0
		for _, p := range occ {
			e1 += hso.At(p, p)
		}

		// Two-electron diagonal
		e2 := 0.0
		for _, p := range occ {
			for _, q := range occ {
				e2 += 0.5 * eri[p][q][p][q]
			}
		}

		H.Set(I, I, e1+e2)
	}

	// Off-diagonal elements
	for I, d := range dets {
		occ := occIndices(d, nspin)

		// Find virtual orbitals
		vir := []int{}
		for a := 0; a < nspin; a++ {
			if (d>>uint(a))&1 == 0 {
				vir = append(vir, a)
			}
		}

		// Singles: i -> a
		for _, i := range occ {
			for _, a := range vir {
				d1 := d ^ (1 << uint(i)) ^ (1 << uint(a))
				J, exists := idx[d1]
				if !exists || J < I {
					continue
				}

				sgn := phaseSingle(d, i, a)
				val := hso.At(a, i)

				for _, p := range occ {
					if p == i {
						continue
					}
					val += eri[a][p][i][p]
				}

				H.Set(I, J, float64(sgn)*val)
			}
		}

		// Doubles: i,j -> a,b
		for ii := 0; ii < len(occ); ii++ {
			i := occ[ii]
			for jj := ii + 1; jj < len(occ); jj++ {
				j := occ[jj]
				for aa := 0; aa < len(vir); aa++ {
					a := vir[aa]
					for bb := aa + 1; bb < len(vir); bb++ {
						b := vir[bb]
						d2 := d ^ (1 << uint(i)) ^ (1 << uint(j)) ^ (1 << uint(a)) ^ (1 << uint(b))
						J, exists := idx[d2]
						if !exists || J < I {
							continue
						}

						sgn := phaseDouble(d, i, j, a, b)
						val := eri[a][b][i][j]
						H.Set(I, J, float64(sgn)*val)
					}
				}
			}
		}
	}

	// Symmetrize
	for i := 0; i < ndet; i++ {
		for j := i + 1; j < ndet; j++ {
			val := H.At(i, j)
			H.Set(j, i, val)
		}
	}

	// Diagonalize
	var eig mat.EigenSym
	ok := eig.Factorize(H, true)
	if !ok {
		panic("FCI diagonalization failed")
	}

	E := make([]float64, ndet)
	eig.Values(E)

	vecs := mat.NewDense(ndet, ndet, nil)
	eig.VectorsTo(vecs)

	return E[0], E, vecs, dets
}