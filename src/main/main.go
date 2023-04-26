package main

import (
	"fmt"
	"log"
	"time"

	ui "github.com/gizak/termui/v3"
	"github.com/gizak/termui/v3/widgets"
	"github.com/src/src/simulation"
)

// position, aceleration, velocity, mass and energy arrays, vector and variables
var pos [][]float64
var acc [][]float64
var vel [][]float64
var mass []float64
var eK float64
var eP float64

// engine param
var screenWidth, screenHeight int = 512, 512

// body generation params
var numberBody int = 4
var maxMass float64 = 0.00005
var maxVel float64 = 3.0
var maxDim float64 = 3.0

// simulation param
var softening float64 = 0.1
var gravity float64 = 1.0
var dt float64 = 0.05
var tCurrent float64 = 0.0
var tStep int = int(tEnd / dt)
var tEnd float64 = 100

func main() {

	// generates initial condition for the simulation, before the main loop
	pos, vel, mass = simulation.GenerateBodies(numberBody, maxDim, maxMass, maxVel)
	acc = simulation.CalcAcc(pos, mass, numberBody, gravity, softening)
	simulation.CalcVel(&vel, acc, numberBody, dt)
	simulation.CalcPos(&pos, vel, numberBody, dt)
	eK, eP = simulation.CalcForces(pos, mass, vel, numberBody, gravity)

	//plotting data
	if err := ui.Init(); err != nil {
		log.Fatal(err)
	}
	defer ui.Close()

	t := widgets.NewParagraph()
	t.Text = "PRESS q TO QUIT DEMO"
	t.SetRect(0, 0, 50, 3)
	t.TextStyle.Fg = ui.ColorWhite
	t.BorderStyle.Fg = ui.ColorCyan

	f := widgets.NewParagraph()
	f.SetRect(0, 3, 50, 8)
	f.TextStyle.Fg = ui.ColorWhite
	f.BorderStyle.Fg = ui.ColorCyan

	p := widgets.NewPlot()
	p.Title = "n-body simulation"
	p.Marker = widgets.MarkerDot
	p.ShowAxes = true

	p.SetRect(0, 10, 50, 30)
	p.Data = pos
	p.AxesColor = ui.ColorWhite
	p.LineColors[0] = ui.ColorWhite
	p.LineColors[1] = ui.ColorWhite
	p.LineColors[2] = ui.ColorWhite
	p.PlotType = widgets.ScatterPlot

	uiEvents := ui.PollEvents()

	draw := func(count int) {
		ui.Render(t, f, p)
	}

	tickerCount := 1
	draw(tickerCount)
	tickerCount++
	ticker := time.NewTicker(time.Millisecond).C

	for {
		f.Text = fmt.Sprintf("ek : %f, ep : %f, eTot : %f, current time : %f\r", eK, eP, eK-eP, tCurrent)
		select {
		case e := <-uiEvents:
			switch e.ID {
			case "q", "<C-c>":
				return
			}
		case <-ticker:
			draw(tickerCount)
			tickerCount++
		}

		simulation.CalcVel(&vel, acc, numberBody, dt)
		simulation.CalcPos(&pos, vel, numberBody, dt)
		acc = simulation.CalcAcc(pos, mass, numberBody, gravity, softening)
		simulation.CalcVel(&vel, acc, numberBody, dt)

		eK, eP = simulation.CalcForces(pos, mass, vel, numberBody, gravity)

		tCurrent += dt

	}

	//p := plot.New()
	//scatter, _ := plotter.NewScatter(plottingq.MakePoints(pos, numberBody))
	// loops in time steps looping until the end time is reached

	//scatter, _ = plotter.NewScatter(plotting.MakePoints(pos, numberBody))
	//p.Add(scatter)

	//fmt.Printf("ek : %f, ep : %f, eTot : %f, current time : %f\r", eK, eP, eK+eP, tCurrent)

	//p.Save(512, 512, "./images/image.png")

}
