(ns stat_prim
	(import Jama.Matrix))

(defn mean [vec] 
	(/ (reduce + vec) (count vec)))

(defn squared [x] (* x x))

(defn disp [vec] 
	(let [m (mean vec)]
		 (/ (reduce + (map #(squared (- %1 m)) vec)) (- (count vec) 1))))

(defn sd [vec] (Math/sqrt (disp vec)))

(defn coef [veca vecb] 
	(let [ma (mean veca)
		  mb (mean vecb)
		  a (reduce + (map * (map #(- %1 ma) veca) (map #(- %1 mb) vecb)))
		  b (* (reduce + (map #(squared (- %1 ma)) veca)) (reduce + (map #(squared (- %1 mb)) vecb)))]
		(/ a (Math/sqrt b))))

	(defn cov [veca vecb] 
		(let [meana (mean veca)
			  meanb (mean vecb)]
		(/ (reduce + (map * (map #(- %1 meana) veca) (map #(- %1 meanb) vecb))) (count veca))))

(defn simple-reg [veca vecb] 
	(let [vecab (map * veca vecb)
		  m1 (- (mean vecab) (* (mean veca) (mean vecb)))
		  m2 (- (mean (map squared veca)) (squared (mean veca)))
		  b (/ m1 m2)
		  a (- (mean vecb) (* b (mean veca)))]
		[a b]))

(defn normal-den-func [x & {:keys [d u] :or {d 1 u 0}}] 
	(/ (Math/exp (- (/ (squared (- x u)) (* 2 d)))) (Math/sqrt (* 2 Math/PI d))))


(defn erf [x] 
	(let [ a1  0.254829592
           a2  -0.284496736
           a3   1.421413741
           a4  -1.453152027 
           a5   1.061405429
           p    0.3275911
		   sign  (if (< x 0) -1 1)
		   xi  (Math/abs x)
		   t  (/ 1 (+ 1 (* p xi)))
		   y  (- 1 (*  (+ a1 (* t (+ a2 (* t (+ a3 (* t (+ a4 (* a5 t)))))))) t (Math/exp (- (squared xi)))))
		]
		(* sign y)))
(defn ercf [x] (- 1 (erf x)))
(defn laplas [x] (* 1/2 (erf (/ x (Math/sqrt 2)))))

(defn normal-dis-func [x & {:keys [d u] :or {d 1 u 0}}] 
	(/ (+ 1 (erf (/ (- x u) (Math/sqrt (* d 2))))) 2))


(defn chi-square-dis [a n] 
	(let [d (if (>= a 0.5) (+ (* 2.0637 (Math/pow (- (Math/log (/ 1 (- 1 a))) 0.16) 0.4274) ) -1.5774)
						   (+ (* -2.0637 (Math/pow (- (Math/log (/ 1 (- 1 a))) 0.16) 0.4274))1.5774))
		  A (* d (Math/sqrt 2))
		  B (* 2/3 (- (squared d) 1))
		  C (* d (/ (- (squared d) 7) (* 9 (Math/sqrt 2))))
		  D (/ (+ (* 6 (squared (squared d))) (* 14 (squared d)) -32) 405)
		  E (* d (/ (+ (* 4 (squared (squared d))) (* 256 (squared d)) -433) (* 4860 (Math/sqrt 2))))]
		(+ n (* A (Math/sqrt n)) B (/ C (Math/sqrt n)) (/ D n) (/ E (* n (Math/sqrt n))))))

(defn chi-square-dis2 [a n] 
	(let [d (if (>= a 0.5) (+ (* 2.0637 (Math/pow (- (Math/log (/ 1 (- 1 a))) 0.16) 0.4274) ) -1.5774)
						   (+ (* -2.0637 (Math/pow (- (Math/log (/ 1 (- 1 a))) 0.16) 0.4274))1.5774))
		  a [1.0000886 0.4713941 0.0001348028 -0.008553069 0.00312558 -0.0008426812 0.00009780499]
		  b [-0.2237368 0.02607083 0.01128186 -0.01153761 0.005169654 0.00253001 -0.001450117]
		  c [-0.001450117 -0.008986007 0.02277679 -0.01323293 -0.006950356 0.001060438 0.001565326]]
		(* n (Math/pow (reduce + (map #(* (Math/pow n (- (/ %1 2))) (Math/pow d %1) (+ %2 (/ %3 n) (/ %4 (squared n)))) (range 0 7) a b c)) 3))))


(defn intervals-split [myvec n] 
	(let [xmin (reduce min myvec)
		  xmax (reduce max myvec)
		  step (inc (quot (- xmax xmin) n))
		  col (vec (doall(repeat n '())))
		  ]
	 (loop 	[v myvec acc col]
		(if	 (empty? v) (with-meta (filter #(not (empty? %)) acc) {:step step})
			(let [x (peek v)
				  pos (if (zero? (mod x step)) (dec (quot x step)) (quot x step)) ]
			(recur (pop v) (assoc acc  (if (>= pos  0) pos 0) (cons x (get acc pos)))))))))

(defn interval-freq [vals] 
	(let [all (reduce concat vals)
		  len (count all)]
	(map #(/ (count %) len) vals)))

(defn p-value [chi2 n] 
	(let [first-stage-range (range 0.01 1 0.1)
		  first-stage-dist (map #(Math/abs (- (chi-square-dis % n) chi2)) first-stage-range)
		  first-stage-i (first (first (sort-by first (map #(vec [%1 %2]) first-stage-range first-stage-dist))))
		  second-stage-range (range first-stage-i (+ 0.1 first-stage-i) 0.01)
		  second-stage-dist (map #(Math/abs (- (chi-square-dis % n) chi2)) second-stage-range)
		  second-stage-i (first (first (sort-by first (map #(vec [%1 %2]) second-stage-range second-stage-dist))))]
		  (nth second-stage-range second-stage-i)))


(defn uniform-test [myvec & {:keys [o] :or {o 0.975}}]  
	(let [intervals (filter #(not (empty? %)) (intervals-split myvec 20))
		  n-freq (/ 1 (count  intervals))
		  t-count (* n-freq (count myvec))
		  chi2 (reduce + (map #(/ (squared (- % t-count)) t-count) (map count intervals)))]
		[ (< chi2 (chi-square-dis o (dec (count intervals)))) chi2 ]))
		
(defn normal-test [myvec &{:keys [o] :or {o 0.975}}]
	(let [intervals (intervals-split myvec 20)
		 n (count myvec)
		 m (mean myvec)
		 d (disp myvec)
		 vecmin (reduce min myvec)
		 vecmax (reduce max myvec)
		 istep (get (meta intervals) :step)
		 interval-borders (for [x (range (count intervals))] [(+ vecmin  (* x istep)) (dec (+ istep vecmin  (* x istep)))])
		 interval-means (map mean interval-borders)
		 t-freqs (map #(- (laplas (float (/ (- (+ % istep) m) d )))  (laplas (float (/ (- % m) d )) )) interval-means) 
		 t-count (map #(* n %) t-freqs)
		 chi2 (reduce + (map #(/ (squared (- %1 %2)) %2) (map count intervals) t-count))]
		  [ (< chi2 (chi-square-dis o (dec (count intervals)))) chi2 ])) 
			
(defn r-parameters [resvec matrix model ]
	(let [ymean (mean resvec)
		 ytr (map #(reduce + (map * % model)) matrix)
		 tss (reduce + (map #(squared (- % ymean)) resvec))
		 rss (reduce + (map #(squared (- %1 %2)) resvec ytr))
         m  (Matrix. (into-array (map double-array matrix)))
         mt (map vec (-> m .transpose (.times m) .inverse .getArray))
		 r-squared (- 1 (/ rss tss))
		 n (dec (count matrix))
		 k (count (first matrix))
         r-adjust (- 1 (* (- 1 r-squared) (/ (- n 1) (- n k))))
         stde (map #(Math/sqrt (* (/ rss (- n k 1))  (nth (nth mt %) %))) (range (count model)))
         t-stat (map #(if (zero? %2) 0 (/ %1 %2)) model stde)]
         {:r-squared r-squared :r-adjusted r-adjust :stde stde :t-stat t-stat})) 			

(defn linear-fit [resvec matrix ] 
	 (let [matrix+const (map #(reverse (conj (reverse %) 1)) (vec (doall matrix)))
		  md (into-array (map double-array matrix+const))
		  rd (double-array resvec)
		  m (Matrix. md)
		  r (Matrix. rd 1)
		  model (flatten (map vec (try (-> m .transpose (.times m) .inverse (.times (.transpose m)) (.times (.transpose r)) .getArray)
	 (catch Exception e (-> m .transpose (.times (-> m (.times (.transpose m)) .inverse)) (.times (.transpose r)) .getArray)))))]

        (with-meta model (r-parameters resvec matrix+const model ))))


	 




	
		   



