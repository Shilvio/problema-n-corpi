package simulation

import (
	"math"
	"math/rand"
)

// generate random bodies all params can be setted in the main function
func GenerateBodies(n int, dim float64, m float64, v float64) ([][]float64, [][]float64, []float64) {

	pos := make([][]float64, 3)
	vel := make([][]float64, 3)

	for i := 0; i < 3; i++ {
		pos[i] = make([]float64, n)
		vel[i] = make([]float64, n)
	}

	mass := make([]float64, n)
	for i := 0; i < n; i++ {

		mass[i] = rand.Float64() * m

		pos[0][i] = rand.NormFloat64() * dim
		pos[1][i] = rand.NormFloat64() * dim
		pos[2][i] = rand.NormFloat64() * dim

		vel[0][i] = rand.NormFloat64() * v
		vel[1][i] = rand.NormFloat64() * v
		vel[2][i] = rand.NormFloat64() * v

	}
	return pos, vel, mass

}

/*
calculates acceleration using Newton's gravitational law, a softening factor is used when calculating distances
so that the function doesn't go to infinite while calculating particles too close together
*/
func CalcAcc(pos [][]float64, mass []float64, n int, g float64, softening float64) [][]float64 {

	acc := make([][]float64, 3)

	for i := 0; i < 3; i++ {
		acc[i] = make([]float64, n)
	}

	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			if i == j {
				continue
			}
			dx := pos[0][i] - pos[0][j]
			dy := pos[1][i] - pos[1][j]
			dz := pos[2][i] - pos[2][j]

			inv_r3 := math.Pow((dx*dx + dy*dy + dz*dz + softening*softening), -1.5)

			acc[0][i] = g * dx * inv_r3 * mass[j]
			acc[1][i] = g * dy * inv_r3 * mass[j]
			acc[2][i] = g * dz * inv_r3 * mass[j]

		}
	}
	return acc
}

// only used to validate the system, calculates potential and kinetic energy
func CalcForces(pos [][]float64, mass []float64, vel [][]float64, n int, g float64) (float64, float64) {

	var eK, eP float64 = 0, 0
	for i := 0; i < n; i++ {

		velVec := math.Sqrt(vel[0][i]*vel[0][i] + vel[1][i]*vel[1][i] + vel[2][i]*vel[2][i])

		for j := 0; j < n; j++ {

			if j == i {
				continue
			}

			dx := pos[0][i] - pos[0][j]
			dy := pos[1][i] - pos[1][j]
			dz := pos[2][i] - pos[2][j]

			inv_r := math.Sqrt(dx*dx + dy*dy + dz*dz)
			inv_r = 1.0 / inv_r
			eK += 0.5 * mass[i] * velVec
			eP += g * -(mass[j] * mass[i]) * inv_r

		}
	}
	return eK, eP
}

// calculate half step velocity
func CalcVel(vel *[][]float64, acc [][]float64, n int, dt float64) {
	v := *vel
	for i := 0; i < n; i++ {
		v[0][i] += acc[0][i] * dt / 2.0
		v[1][i] += acc[1][i] * dt / 2.0
		v[2][i] += acc[2][i] * dt / 2.0
	}
}

// calculate full step position
func CalcPos(pos *[][]float64, vel [][]float64, n int, dt float64) {
	p := *pos
	for i := 0; i < n; i++ {
		p[0][i] += vel[0][i] * dt
		p[1][i] += vel[1][i] * dt
		p[2][i] += vel[2][i] * dt
	}

}
