package main

import (
	"GoProject/BasicOneD"
	"GoProject/OurServer"
	"log"
)

func main() {
	OurServer.WebServer()
	err := BasicOneD.OneDimPoissonSolver()
	if err != nil {
		log.Fatal(err)
	}

	err = BasicOneD.BarrierPotential()
	if err != nil {
		log.Fatal(err)
	}
}
