# CLojure statistical primitives

Usable, but not useful yet.


## Example

$ (mean [1 2 3 4 5]) => 3
$ (uniform-test [1 1 1 1 1 0 0 0 0 0]) =>  [true 0N]
$ (normal-test [1 2 2 3 3 3 4 4 5]) => [true 0.8779279346018437]
$ (def c (linear-fit [1 2 3 4 5] (map #(vector (* 0.5 %)) [1 2 3 4 5])))
$ c => (1.9999999999999996 1.9984014443252818E-15)
$ (meta c) => {:r-squared 1.0, :r-adjusted 1.0, :stde (1.8630587706705473E-15 3.0895334523474967E-15), :t-stat (1.0735034404095394E15 0.6468295213980776)}







