package traceskeleton

import (
	"fmt"
	"sync"
)

func abs(x int) int {
	if x < 0 {
		return -x
	}
	return x
}

//================================
// ENUMS
//================================

const (
	HORIZONTAL = 1
	VERTICAL   = 2
)

//================================
// DATASTRUCTURES
//================================

type Point struct {
	X    int
	Y    int
	Next *Point
}
type Polyline struct {
	Head *Point
	Tail *Point
	Prev *Polyline
	Next *Polyline
	Size int
}

//================================
// DATASTRUCTURE IMPLEMENTATION
//================================

func newPolyline() *Polyline {
	var q *Polyline = new(Polyline)
	q.Head = nil
	q.Tail = nil
	q.Prev = nil
	q.Next = nil
	q.Size = 0
	return q
}
func PrintPolyline(q *Polyline) {
	if q == nil {
		return
	}
	var jt *Point = q.Head
	for jt != nil {
		fmt.Printf("%d,%d ", jt.X, jt.Y)
		jt = jt.Next
	}
	fmt.Printf("\n")
}
func PrintPolylines(q *Polyline) {
	if q == nil {
		return
	}
	var it *Polyline = q
	for it != nil {
		var jt *Point = it.Head
		for jt != nil {
			fmt.Printf("%d,%d ", jt.X, jt.Y)
			jt = jt.Next
		}
		fmt.Printf("\n")
		it = it.Next
	}
}
func reversePolyline(q *Polyline) {
	if q == nil || q.Size < 2 {
		return
	}
	q.Tail.Next = q.Head
	var it0 *Point = q.Head
	var it1 *Point = it0.Next
	var it2 *Point = it1.Next
	for i := 0; i < q.Size-1; i++ {
		it1.Next = it0
		it0 = it1
		it1 = it2
		it2 = it2.Next
	}
	var qHead *Point = q.Head
	q.Head = q.Tail
	q.Tail = qHead
	q.Tail.Next = nil
}

func catTailPolyline(q0 *Polyline, q1 *Polyline) {
	if q1 == nil {
		return
	}
	if q0 == nil {
		*q0 = *newPolyline()
	}
	if q0.Head == nil {
		q0.Head = q1.Head
		q0.Tail = q1.Tail
		return
	}
	q0.Tail.Next = q1.Head
	q0.Tail = q1.Tail
	q0.Size += q1.Size
	q0.Tail.Next = nil
}

func catHeadPolyline(q0 *Polyline, q1 *Polyline) {
	if q1 == nil {
		return
	}
	if q0 == nil {
		*q0 = *newPolyline()
	}
	if q1.Head == nil {
		return
	}
	if q0.Head == nil {
		q0.Head = q1.Head
		q0.Tail = q1.Tail
		return
	}
	q1.Tail.Next = q0.Head
	q0.Head = q1.Head
	q0.Size += q1.Size
	q0.Tail.Next = nil
}

func addPointToPolyline(q *Polyline, x int, y int) {
	var p *Point = new(Point)
	p.X = x
	p.Y = y
	p.Next = nil
	if q.Head == nil {
		q.Head = p
		q.Tail = p
	} else {
		q.Tail.Next = p
		q.Tail = p
	}
	q.Size++
}

func prependPolyline(q0 *Polyline, q1 *Polyline) *Polyline {
	if q0 == nil {
		return q1
	}
	q1.Next = q0
	q0.Prev = q1
	return q1
}

//================================
// RASTER SKELETONIZATION
//================================
// Binary image thinning (skeletonization) in-place.
// Implements Zhang-Suen algorithm.
// http://agcggs680.pbworks.com/f/Zhan-SuenAlgorithm.pdf

func thinningZSIteration(im []uint8, w int, h int, iter int) bool {
	diff := false
	for i := 1; i < h-1; i++ {
		for j := 1; j < w-1; j++ {
			p2 := im[(i-1)*w+j] & 1
			p3 := im[(i-1)*w+j+1] & 1
			p4 := im[(i)*w+j+1] & 1
			p5 := im[(i+1)*w+j+1] & 1
			p6 := im[(i+1)*w+j] & 1
			p7 := im[(i+1)*w+j-1] & 1
			p8 := im[(i)*w+j-1] & 1
			p9 := im[(i-1)*w+j-1] & 1
			A := 0
			if p2 == 0 && p3 == 1 {
				A++
			}
			if p3 == 0 && p4 == 1 {
				A++
			}
			if p4 == 0 && p5 == 1 {
				A++
			}
			if p5 == 0 && p6 == 1 {
				A++
			}
			if p6 == 0 && p7 == 1 {
				A++
			}
			if p7 == 0 && p8 == 1 {
				A++
			}
			if p8 == 0 && p9 == 1 {
				A++
			}
			if p9 == 0 && p2 == 1 {
				A++
			}
			B := p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9
			var m1 uint8 = 0
			var m2 uint8 = 0
			if iter == 0 {
				m1 = p2 * p4 * p6
				m2 = p4 * p6 * p8
			} else {
				m1 = p2 * p4 * p8
				m2 = p2 * p6 * p8
			}
			if A == 1 && (B >= 2 && B <= 6) && m1 == 0 && m2 == 0 {
				im[i*w+j] |= 2
			}
		}
	}
	for i := 0; i < h*w; i++ {
		var marker bool = (im[i] >> 1) != 0
		var old uint8 = (im[i] & 1)
		if old != 0 && (!marker) {
			im[i] = 1
		} else {
			im[i] = 0
		}
		if (!diff) && (im[i] != old) {
			diff = true
		}
	}
	return diff
}

func ThinningZS(im []uint8, w int, h int) {
	var diff bool = true
	for {
		diff = diff && thinningZSIteration(im, w, h, 0)
		diff = diff && thinningZSIteration(im, w, h, 1)
		if !diff {
			break
		}
	}
}

//================================
// MAIN ALGORITHM
//================================

// check if a region has any white pixel
func notEmpty(im []uint8, W int, H int, x int, y int, w int, h int) bool {
	for i := y; i < y+h; i++ {
		for j := x; j < x+w; j++ {
			if im[i*W+j] != 0 {
				return true
			}
		}
	}
	return false
}

/**merge ith fragment of second chunk to first chunk
 * @param c0   fragments from  first  chunk
 * @param c1i  ith fragment of second chunk
 * @param sx   (x or y) coordinate of the seam
 * @param isv  is vertical, not horizontal?
 * @param mode 2-bit flag,
 *             MSB = is matching the left (not right) end of the fragment from first  chunk
 *             LSB = is matching the right (not left) end of the fragment from second chunk
 * @return     matching successful?
 */
func mergeImpl(c0 *Polyline, c1i *Polyline, sx int, isv bool, mode int) bool {
	var b0 bool = (mode >> 1 & 1) > 0 // match c0 left
	var b1 bool = (mode >> 0 & 1) > 0 // match c1 left
	var c0j *Polyline = nil
	var md int = 4
	var p1 *Point = c1i.Tail
	if b1 {
		p1 = c1i.Head
	}
	var xx int = p1.X
	if isv {
		xx = p1.Y
	}
	if abs(xx-sx) > 0 { // not on the seam, skip
		return false
	}
	// find the best match
	var it *Polyline = c0
	for it != nil {
		var p0 *Point = it.Tail
		if b0 {
			p0 = it.Head
		}
		var aa int = p0.X
		if isv {
			aa = p0.Y
		}
		if abs(aa-sx) > 1 {
			it = it.Next
			continue
		}
		var bb int = p0.Y
		if isv {
			bb = p0.X
		}
		var cc int = p1.Y
		if isv {
			cc = p1.X
		}
		var d int = abs(bb - cc)
		if d < md {
			c0j = it
			md = d
		}
		it = it.Next
	}
	if c0j != nil { // best match is good enough, merge them
		if b0 && b1 {
			reversePolyline(c1i)
			catHeadPolyline(c0j, c1i)
		} else if !b0 && b1 {
			catTailPolyline(c0j, c1i)
		} else if b0 && !b1 {
			catHeadPolyline(c0j, c1i)
		} else {
			reversePolyline(c1i)
			catTailPolyline(c0j, c1i)
		}
		return true
	}
	return false
}

/**merge fragments from two chunks
 * @param c0   fragments from first  chunk
 * @param c1   fragments from second chunk
 * @param sx   (x or y) coordinate of the seam
 * @param dr   merge direction, HORIZONTAL or VERTICAL?
 */
func mergeFrags(c0 *Polyline, c1 *Polyline, sx int, dr int) *Polyline {
	if c0 == nil {
		return c1
	}
	if c1 == nil {
		return c0
	}
	var it *Polyline = c1
	for it != nil {
		var tmp *Polyline = it.Next
		if dr == HORIZONTAL {
			if mergeImpl(c0, it, sx, false, 1) {
				goto rem
			}
			if mergeImpl(c0, it, sx, false, 3) {
				goto rem
			}
			if mergeImpl(c0, it, sx, false, 0) {
				goto rem
			}
			if mergeImpl(c0, it, sx, false, 2) {
				goto rem
			}
		} else {
			if mergeImpl(c0, it, sx, true, 1) {
				goto rem
			}
			if mergeImpl(c0, it, sx, true, 3) {
				goto rem
			}
			if mergeImpl(c0, it, sx, true, 0) {
				goto rem
			}
			if mergeImpl(c0, it, sx, true, 2) {
				goto rem
			}
		}
		goto next
	rem:
		if it.Prev == nil {
			c1 = it.Next
			if it.Next != nil {
				it.Next.Prev = nil
			}
		} else {
			it.Prev.Next = it.Next
			if it.Next != nil {
				it.Next.Prev = it.Prev
			}
		}
	next:
		it = tmp
	}
	it = c1
	for it != nil {
		var tmp *Polyline = it.Next
		it.Prev = nil
		it.Next = nil
		c0 = prependPolyline(c0, it)
		it = tmp
	}
	return c0
}

/**recursive bottom: turn chunk into polyline fragments;
 * look around on 4 edges of the chunk, and identify the "outgoing" pixels;
 * add segments connecting these pixels to center of chunk;
 * apply heuristics to adjust center of chunk
 *
 * @param x    left of   chunk
 * @param y    top of    chunk
 * @param w    width of  chunk
 * @param h    height of chunk
 * @return     the polyline fragments
 */
func chunkToFrags(im []uint8, W int, H int, x int, y int, w int, h int) *Polyline {
	var frags *Polyline = nil
	var fsize int = 0
	var on bool = false // to deal with strokes thicker than 1px
	var li int = -1
	var lj int = -1
	// walk around the edge clockwise
	for k := 0; k < h+h+w+w-4; k++ {
		var i, j int
		if k < w {
			i = y + 0
			j = x + k
		} else if k < w+h-1 {
			i = y + k - w + 1
			j = x + w - 1
		} else if k < w+h+w-2 {
			i = y + h - 1
			j = x + w - (k - w - h + 3)
		} else {
			i = y + h - (k - w - h - w + 4)
			j = x + 0
		}
		if im[i*W+j] != 0 {
			if !on {
				on = true
				var f *Polyline = newPolyline()
				addPointToPolyline(f, j, i)
				addPointToPolyline(f, x+w/2, y+h/2)
				frags = prependPolyline(frags, f)
				fsize++
			}
		} else {
			if on { // right side of stroke, average to get center of stroke
				frags.Head.X = (frags.Head.X + lj) / 2
				frags.Head.Y = (frags.Head.Y + li) / 2
				on = false
			}
		}
		li = i
		lj = j
	}
	if fsize == 2 {
		var f *Polyline = newPolyline()
		addPointToPolyline(f, frags.Head.X, frags.Head.Y)
		addPointToPolyline(f, frags.Next.Head.X, frags.Next.Head.Y)
		frags = f
	} else if fsize > 2 {
		var ms int = 0
		var mi int = -1
		var mj int = -1
		// use convolution to find brightest blob
		for i := y + 1; i < y+h-1; i++ {
			for j := x + 1; j < x+w-1; j++ {
				var s int = int((im[i*W-W+j-1]) + (im[i*W-W+j]) + (im[i*W-W+j-1+1]) +
					(im[i*W+j-1]) + (im[i*W+j]) + (im[i*W+j+1]) +
					(im[i*W+W+j-1]) + (im[i*W+W+j]) + (im[i*W+W+j+1]))
				if s > ms {
					mi = i
					mj = j
					ms = s
				} else if s == ms && abs(j-(x+w/2))+abs(i-(y+h/2)) < abs(mj-(x+w/2))+abs(mi-(y+h/2)) {
					mi = i
					mj = j
					ms = s
				}
			}
		}
		if mi != -1 {
			var it *Polyline = frags
			for it != nil {
				it.Tail.X = mj
				it.Tail.Y = mi
				it = it.Next
			}
		}
	}
	return frags
}

/**Trace skeleton from thinning result.
 * Algorithm:
 * 1. if chunk size is small enough, reach recursive bottom and turn it into segments
 * 2. attempt to split the chunk into 2 smaller chunks, either horizontall or vertically;
 *    find the best "seam" to carve along, and avoid possible degenerate cases
 * 3. recurse on each chunk, and merge their segments
 *
 * @param x       left of   chunk
 * @param y       top of    chunk
 * @param w       width of  chunk
 * @param h       height of chunk
 * @param iter    current iteration
 * @return        an array of polylines
 */
func TraceSkeleton(im []uint8, W int, H int, x int, y int, w int, h int, chunkSize int, maxIter int) *Polyline {
	var frags *Polyline = nil

	if maxIter <= 0 { // gameover
		return frags
	}
	if w <= chunkSize && h <= chunkSize {
		frags = chunkToFrags(im, W, H, x, y, w, h)
		return frags
	}
	var ms int = W + H // number of white pixels on the seam, less the better
	var mi int = -1
	var mj int = -1
	if h > chunkSize {
		for i := y + 3; i < y+h-3; i++ {
			if im[i*W+x] != 0 || im[(i-1)*W+x] != 0 || im[i*W+x+w-1] != 0 || im[(i-1)*W+x+w-1] != 0 {
				continue
			}
			var s int = 0
			for j := x; j < x+w; j++ {
				s += int(im[i*W+j])
				s += int(im[(i-1)*W+j])
			}
			if s < ms {
				ms = s
				mi = i
			} else if s == ms && abs(i-(y+h/2)) < abs(mi-(y+h/2)) {
				// if there is a draw (very common), we want the seam to be near the middle
				// to balance the divide and conquer tree
				ms = s
				mi = i
			}
		}
	}
	if w > chunkSize {
		for j := x + 3; j < x+w-3; j++ {
			if im[W*y+j] != 0 || im[W*(y+h)-W+j] != 0 || im[W*y+j-1] != 0 || im[W*(y+h)-W+j-1] != 0 {
				continue
			}
			var s int = 0
			for i := y; i < y+h; i++ {
				s += int(im[i*W+j])
				s += int(im[i*W+j-1])
			}
			if s < ms {
				ms = s
				mi = -1 // horizontal seam is defeated
				mj = j
			} else if s == ms && abs(j-(x+w/2)) < abs(mj-(x+w/2)) {
				ms = s
				mi = -1
				mj = j
			}
		}
	}
	var L0 int = -1
	var L1, L2, L3 int
	var R0 int = -1
	var R1, R2, R3 int
	var dr int = 0
	var sx int
	if h > chunkSize && mi != -1 { // split top and bottom
		L0 = x
		L1 = y
		L2 = w
		L3 = mi - y
		R0 = x
		R1 = mi
		R2 = w
		R3 = y + h - mi
		dr = VERTICAL
		sx = mi
	} else if w > chunkSize && mj != -1 { // split left and right
		L0 = x
		L1 = y
		L2 = mj - x
		L3 = h
		R0 = mj
		R1 = y
		R2 = x + w - mj
		R3 = h
		dr = HORIZONTAL
		sx = mj
	}
	var aL bool = false
	var aR bool = false
	if dr != 0 && notEmpty(im, W, H, L0, L1, L2, L3) { // if there are no white pixels, don't waste time
		aL = true
	}
	if dr != 0 && notEmpty(im, W, H, R0, R1, R2, R3) {
		aR = true
	}
	if aL && aR {
		messages := make(chan int)
		var wg sync.WaitGroup
		wg.Add(2)
		var pL, pR *Polyline
		go func() {
			defer wg.Done()
			pL = TraceSkeleton(im, W, H, L0, L1, L2, L3, chunkSize, maxIter-1)
			messages <- 1
		}()
		go func() {
			defer wg.Done()
			pR = TraceSkeleton(im, W, H, R0, R1, R2, R3, chunkSize, maxIter-1)
			messages <- 2

		}()
		go func() {
			wg.Wait()
			close(messages)
		}()
		<-messages
		<-messages
		frags = mergeFrags(pL, pR, sx, dr)

		// no goroutines
		// frags = mergeFrags(
		//     traceSkeleton(im,W,H,L0,L1,L2,L3,chunkSize,maxIter-1),
		//     traceSkeleton(im,W,H,R0,R1,R2,R3,chunkSize,maxIter-1),
		// sx,dr)
	} else if aL {
		frags = TraceSkeleton(im, W, H, L0, L1, L2, L3, chunkSize, maxIter-1)
	} else if aR {
		frags = TraceSkeleton(im, W, H, R0, R1, R2, R3, chunkSize, maxIter-1)
	}
	if mi == -1 && mj == -1 { // splitting failed! do the recursive bottom instead
		frags = chunkToFrags(im, W, H, x, y, w, h)
	}
	return frags
}

func PolylinesToSvg(q *Polyline, w int, h int) string {
	var svg string = fmt.Sprintf("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%d\" height=\"%d\" fill=\"none\" stroke=\"black\" stroke-width=\"1\">", w, h)
	if q == nil {
		return svg + "</svg>"
	}
	var it *Polyline = q
	for it != nil {
		var jt *Point = it.Head
		svg += "<path d=\"M"
		for jt != nil {
			svg += fmt.Sprintf("%d,%d", jt.X, jt.Y)
			jt = jt.Next
			if jt != nil {
				svg += " L"
			}
		}
		svg += "\"/>"
		it = it.Next
	}
	return svg + "</svg>"
}
