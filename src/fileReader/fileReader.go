package fileReader

import (
	"bufio"
	"fmt"
	"log"
	"os"
)

type Xy struct {
	Id int
	X  float64
	Y  float64
}

func readData(path string) ([]Xy, error) {

	file, err := os.Open(path)
	if err != nil {
		log.Fatalf("could not open file %v", err)
	}
	defer file.Close()

	var value []Xy
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		var id int
		var x, y float64

		i, err := fmt.Sscanf(scanner.Text(), "%d,%f,%f", &id, &x, &y)
		if i == 0 {
			continue
		}
		if err != nil {
			log.Println("bad data point ", scanner.Text(), " : %v", err)
		}

		value = append(value, Xy{id, x, y})
	}
	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("could not read file %v", err)
	}
	return value, nil
}

func ParseFile() []Xy {
	value, err := readData("../../csv/particle.csv")
	if err != nil {
		log.Fatal("could not read file")
	}
	_ = value

	return value
}
