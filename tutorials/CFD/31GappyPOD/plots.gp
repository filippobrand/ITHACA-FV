# Output configuration (saves as PNG; change to 'set terminal qt' for an interactive window)
set terminal pngcairo size 1000,600 font "Helvetica,11"
set output 'l2_error_comparison.png'

# Titles and Axis Labels
set title "Gappy reconstruction error comparison" font ",13"
set xlabel "Sample Index"
set ylabel "L_2 Relative error" font ",12"

# Grid and Visual Styling
set grid
set key top right spacing 1.2
set tics nomirror

# Define line styles for clear distinction
set style line 1 linecolor rgb "#E69F00" lw 2 pt 7 ps 0.6  # Orange
set style line 2 linecolor rgb "#56B4E9" lw 2 pt 5 ps 0.6  # Sky Blue
set style line 3 linecolor rgb "#009E73" lw 2 pt 9 ps 0.6  # Bluish Green

# Define base path
dir = "./ITHACAoutput/l2_error/"

# Plotting using zero-based line index ($0) as the X-axis
plot dir."l2_error_3modes_30samples_mat.txt" using 0:1 with linespoints ls 1 title "3 Modes", \
     dir."l2_error_4modes_30samples_mat.txt" using 0:1 with linespoints ls 2 title "4 Modes", \
     dir."l2_error_5modes_30samples_mat.txt" using 0:1 with linespoints ls 3 title "5 Modes"