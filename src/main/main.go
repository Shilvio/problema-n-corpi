package main

import (
	"fmt"

	"github.com/src/src/simulation"
)

// position, aceleration, velocity, mass and energy arrays, vector and variables
var pos [][]float64
var acc [][]float64
var vel [][]float64
var mass []float64
var eK float64
var eP float64

// body generation params
var numberBody int = 2
var maxMass float64 = 0.2
var maxVel float64 = 3.0
var maxDim float64 = 3.0

// simulation param
var softening float64 = 0.1
var gravity float64 = 1.0
var dt float64 = 0.01
var tCurrent float64 = 0
var tStep int = int(tEnd / dt)
var tEnd float64 = 100

func main() {

	pos, vel, mass = simulation.GenerateBodies(numberBody, maxDim, maxMass, maxVel)
	acc = simulation.CalcAcc(pos, mass, numberBody, gravity, softening)

	simulation.CalcVel(&vel, acc, numberBody, dt)
	simulation.CalcPos(&pos, vel, numberBody, dt)
	eK, eP = simulation.CalcForces(pos, mass, vel, numberBody, gravity)

	for i := 0; i < tStep; i++ {

		simulation.CalcVel(&vel, acc, numberBody, dt)
		simulation.CalcPos(&pos, vel, numberBody, dt)
		acc = simulation.CalcAcc(pos, mass, numberBody, gravity, softening)
		simulation.CalcVel(&vel, acc, numberBody, dt)

		eK, eP = simulation.CalcForces(pos, mass, vel, numberBody, gravity)

		tCurrent += dt

		fmt.Printf("ek : %f, ep : %f, eTot : %f, current time : %f\r", eK, eP, eK+eP, tCurrent)
	}

}
