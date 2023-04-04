/* package main

import (
	"visualizer/demo/fileReader"
	graphbuilder "visualizer/demo/graphBuilder"
)

type xy = fileReader.Xy

func main() {

	var positions []xy
	positions = fileReader.ParseFile()
	//for _, a := range positions {
	//	fmt.Println(a.Id, a.X, a.Y)
	//}
	graphbuilder.GenerateGraph(positions)
}
*/

/*



 */

package main

import (
	"C"
)
import "fmt"

func main() {

	fmt.Println("calling c barnes-hut")

}
