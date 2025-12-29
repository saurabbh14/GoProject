package main

import (
	"GoProject/BasicOneD"
	"log"
)

func main() {
	err := BasicOneD.BarrierPotential()
	if err != nil {
		log.Fatal(err)
	}
}
