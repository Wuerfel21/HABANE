# HABANE - **Ha**lf-**ba**nd **N**oise **E**nvelope

This is really just a proof of concept for an idea that has been on my mind for a while.

Essentially, this combines a multiply-free 4-bit ADPCM codec at 24kHz sample rate with a 1.5 kHz 4-bit noise envelope that the decoder uses to approximately restore essential HF content (mostly percussion instruments). 

The goal is to encode reasonable(tm)-quality music at around 240 kbps (= one 512 byte block per 60Hz video frame) while using very little compute and no scratch memory in the decoder. HABANE encodes 434 samples @ 24kHz (27 units of compressed 16 samples + 2 raw 16 bit samples to allow independent decoding of blocks) per block, resulting in a rate of ~55 blocks per second, equalling 226 kbps (or 221 kbps when 1k = 1024, as it should be). SUCCESS! Note that each block contains some currently unused fields intended for looping and conditional branching.

