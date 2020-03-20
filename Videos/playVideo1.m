function h = playVideo1(video_name)

hmfr = video.MultimediaFileReader(video_name);
 hvp = video.VideoPlayer;
 
 while ~isDone(hmfr)
  frame = step(hmfr);
  h = step(hvp, frame);
  h;
 end

 release(hmfr);
 release(hvp);