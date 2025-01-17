#!/bin/bash

# Parameters
CHR="chrX"
CHR_SIZE=124992030  # Replace with the actual size of chrX
CHUNK_SIZE=5000000  # 5 Mb
OVERLAP=100000      # 100 kb

# Output file
OUTPUT="Dog10K.chunks.chrX.txt"

# Initialize chunk start and end
START=1
END=$CHUNK_SIZE

# Generate chunks
> $OUTPUT
while [ $START -lt $CHR_SIZE ]; do
    if [ $END -gt $CHR_SIZE ]; then
        END=$CHR_SIZE
    fi
    echo -e "$CHR $CHR:$START-$END $START $END" >> $OUTPUT
    START=$((START + CHUNK_SIZE - OVERLAP))
    END=$((START + CHUNK_SIZE - 1))
done

echo "Chunks for $CHR written to $OUTPUT"

