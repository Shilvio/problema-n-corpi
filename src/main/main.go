package main

import (
	"fmt"
	"log"
	"os/exec"
)

func main() {

	fmt.Println("calling c barnes-hut")
	out, err := exec.Command("ls", "-l", ".").Output()
	if err != nil {
		log.Fatal("command reported an error ", err)

	}
	fmt.Println("outpout is: %s", out)
}
