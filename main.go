package main

import (
	"GoProject/BasicOneD"
	"GoProject/OurServer"
	"log"
)

func main() {
	OurServer.BasicWebServer()
	err := BasicOneD.BarrierPotential()
	if err != nil {
		log.Fatal(err)
	}
}
