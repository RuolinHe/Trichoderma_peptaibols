function slide_window = Seq2window_HRL(seq,window,slide)
%Seq2window_HRL 为寻找重组区域准备滑窗
if length(seq)<window
    error('sequence length is shorter than window size')
else
    slide_window = cell(ceil((length(seq)-window)/slide)+1,1);
    slide_window{1}=seq(1:window);
    for k = 1:length(slide_window)-2
        slide_window{k+1}=seq(1+k*slide:window+k*slide);
    end
    slide_window{end}=seq(window+k*slide+1:end);
    if length(slide_window{end})<window-slide%最后一个区域太短的话丢掉
        slide_window(end)=[];
    end
end
end

