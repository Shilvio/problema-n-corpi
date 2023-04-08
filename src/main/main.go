package main

import (
	"bufio"
	"fmt"
	"os"
	"strconv"
	"strings"

	"log"
	"os/exec"

	"github.com/hajimehoshi/ebiten/v2"
	"github.com/hajimehoshi/ebiten/v2/ebitenutil"
)

// particle dot
type particle struct {
	id   float64
	xPos float64
	yPos float64
	img  *ebiten.Image
}

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
		dataTmp = append(dataTmp, particle{id, xPos, yPos, ebiten.NewImage(1, 1)})

	}

	return dataTmp
}

// screen options
const screenWidth = 512
const screenHeight = 512

// algorithm calling and particle generation functions
func generateParticlesCall() {
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
	return nil
}

func (g *Game) Draw(screen *ebiten.Image) {
	ebitenutil.DebugPrint(screen, fmt.Sprintf("FPS: %0.2f", ebiten.ActualFPS()))
}

func (g *Game) Layout(outsideWidth, outsideHeight int) (screenWidth, screenHeight int) {
	return 512, 512
}

func main() {
	/*ebiten.SetWindowSize(screenWidth, screenHeight)

	if err := ebiten.RunGame(&Game{}); err != nil {
		log.Fatal(err)
	}*/

	data := loadParticle()
	fmt.Println(data)
}
