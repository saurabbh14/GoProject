package main

import (
	"GoProject/BasicOneD"
	"GoProject/BasicTest"
	"fmt"
	"log"
)

func main() {
	testHarmonicOscillator()
	err := BasicOneD.BarrierPotential()
	if err != nil {
		log.Fatal(err)
	}
}

func testHarmonicOscillator() {
	// Problem parameters
	M := 20
	Ne := 2
	Ngrid := 401
	L := 12.0
	aSoft := 1.0
	strength := 1.0

	// External potential: harmonic oscillator
	vext := func(x float64) float64 {
		return 0.5 * x * x
	}

	fmt.Println("Building integrals...")
	_, h, g := BasicTest.BuildIntegrals(M, Ngrid, L, aSoft, strength, vext)

	fmt.Println("\nRunning RHF...")
	EHF, _, _, _, _ := BasicTest.RHF(h, g, Ne, 100, 1e-10, 0.25, true)
	fmt.Printf("\nRHF energy: %.12f\n", EHF)

	fmt.Println("\nRunning FCI...")
	E0, _, _, _ := BasicTest.FCI(h, g, Ne, true)
	fmt.Printf("FCI ground energy: %.12f\n", E0)
	fmt.Printf("Correlation (FCI - HF): %.12e\n", E0-EHF)
}
