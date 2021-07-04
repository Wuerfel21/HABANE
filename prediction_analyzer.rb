

bestcount = Array.new(512){0}

if (ARGV[0] =~ /\.txt$/)
    puts "parsing file.."
    File.read(ARGV[0]).lines.each do |line|
        line =~ /([0-9A-F]{3}) : (\d+)/
        bestcount[$1.to_i(16)] = $2.to_i
    end
else
    stream = File.binread(ARGV[0]).unpack("s<*")
    stream_l = Array.new
    stream_r = Array.new
    stream.each_slice(2){|(l,r)|stream_l<<l;stream_r<<r}
    stream = stream_l+stream_r

    (2...stream.length).each do |i|

        best_diff = 0xFFFFFFFFFFFFFFF
        best_predicitor = -1

        512.times do |predictor|
            h1sar = (predictor)&15
            h2sar = (predictor>>4)&15
            inverth2 = predictor[8]>0

            h1 = stream[i-1]
            h2 = stream[i-2]
            prediction = h1>>h1sar + (inverth2 ? -h2 : h2)>>h2sar

            diff = (stream[i] - prediction).abs
            if diff < best_diff
                best_diff = diff
                best_predicitor = predictor
            end

        end
        
        bestcount[best_predicitor]+=1
    end
end

#p bestcount
puts "64 best predictors:"
(bestcount.map.with_index.sort.reverse)[0...64].each do |count,predictor|
    h1sar = (predictor)&15
    h2sar = (predictor>>4)&15
    inverth2 = predictor[8]>0
    puts "h1>>%2d, %sh2>>%2d : %d bests" % [h1sar,inverth2 ? ?- : ?+,h2sar,count]
end

