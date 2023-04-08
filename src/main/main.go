package main

import (
	"fmt"
	"log"
	"os/exec"
)

func main() {

	fmt.Println("calling c barnes-hut")
	_, err := exec.Command("gcc", "../barnes-hut-algorithm/Barnes-hut-Bounding-box.c", "-o", "Barnes-hut-Bounding-box").Output()
	if err != nil {
		log.Fatal("command reported an error ", err)

	} else {
		fmt.Println("Program executed")
	}
}
