set -x

render_video () {
    mkdir temp
    cp output/*.png temp/
    convert temp/*.png -delay 10 -morph 10 temp/%05d.jpg
    yes | ffmpeg -r 60 -i temp/%05d.jpg -r 20 output.mp4
    rm -R temp
}

mkdir output
rm output/*
cargo run && render_video
