package main

import (
	"bufio"
	"image/color"

	"fmt"
	"log"
	"math/rand"
	"os"
	"os/exec"
	"strconv"
	"strings"

	"github.com/hajimehoshi/ebiten/v2"
	"github.com/hajimehoshi/ebiten/v2/ebitenutil"
)

// particle structure
type particle struct {
	id   float64
	xPos float64
	yPos float64
	img  *ebiten.Image
}

// global particle's array
var data []particle

func loadParticle() []particle {
	file, err := os.Open("../particles-data/particle.txt")
	if err != nil {
		log.Fatal("File error: ", err)
	}
	dataTmp := []particle{}
	defer file.Close()
	scanner := bufio.NewScanner(file)

	//skip first line
	var skip = false
	for scanner.Scan() {
		if skip == false {
			skip = true
			continue
		}

		//loads variable in particle struct array
		tmp := strings.Split(scanner.Text(), " ")
		id, _ := strconv.ParseFloat(tmp[2], 64)
		xPos, _ := strconv.ParseFloat(tmp[0], 64)
		yPos, _ := strconv.ParseFloat(tmp[1], 64)
		img := ebiten.NewImage(2, 2)
		img.Fill(color.White)
		dataTmp = append(dataTmp, particle{id, xPos, yPos, img})

	}

	return dataTmp
}

// screen options
var screenWidth = 1024
var screenHeight = 1024

// algorithm calling and particle generation functions
func algorithmCall() {
	fmt.Println("calling c barnes-hut")
	_, err := exec.Command("../generate-particles/particleRand").Output()
	if err != nil {
		log.Fatal("command reported an error ", err)

	} else {
		fmt.Println("Program executed")
	}
}

func barnesHutCall() {

	fmt.Println("calling c barnes-hut")
	_, err := exec.Command("../barnes-hut-algorithm/Barnes-hut-Bounding-box").Output()
	if err != nil {
		log.Fatal("command reported an error ", err)

	}

}

// graphics function

// set particle color
func init() {
}

type Game struct{}

func (g *Game) Update() error {

	for _, elem := range data {
		fmt.Println(elem)
	}

	return nil
}

func (g *Game) Draw(screen *ebiten.Image) {
	ebitenutil.DebugPrint(screen, fmt.Sprintf("FPS: %0.2f", ebiten.ActualFPS()))
	options := &ebiten.DrawImageOptions{}

	for _, elem := range data {

		img := elem.img
		options.GeoM.Translate(float64(rand.Intn(512)), float64(rand.Intn(512)))
		screen.DrawImage(img, options)
	}

}

func (g *Game) Layout(outsideWidth, outsideHeight int) (screenWidth, screenHeight int) {
	return screenWidth, screenWidth
}

func main() {
	data := loadParticle()
	if data == nil {
		log.Fatal("invalid data")
		os.Exit(1)
	}
	ebiten.SetWindowSize(screenWidth, screenHeight)
	if err := ebiten.RunGame(&Game{}); err != nil {
		log.Fatal(err)
	}

	//data := loadParticle()
	//if data == nil {
	//log.Fatal("invalid data")
	//}

}
