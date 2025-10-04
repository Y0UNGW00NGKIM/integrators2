outdir=results
mkdir -p "$outdir"

Ns="64 128 256 512 1024 2048 4096 8192 16384 32768 65536 131072 262144 524288 1048576 2097152 4194304 8388608 16777216"

for D in 3 5 10; do
  csv="$outdir/results_d${D}.csv"
  echo "sqrtN,volume,se" > "$csv"
  for N in $Ns; do
    O=$(python3 ndcrescent.py $D $N 1.0 1.0 0.5)
    V=$(echo "$O" | awk '/^volume:/ {print $2}')
    S=$(echo "$O" | awk '/^stat uncertainty:/ {print $3}')
    SQ=$(python3 - << EOF
import math
print(f"{math.sqrt($N):.10f}")
EOF
)
    echo "$SQ,$V,$S" >> "$csv"
  done
done
