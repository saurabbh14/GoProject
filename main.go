package main

import (
	"GoProject/BasicOneD"
	"GoProject/OurServer"
	"log"
)

func main() {
	OurServer.WebServer()
	err := BasicOneD.BarrierPotential()
	if err != nil {
		log.Fatal(err)
	}
}
