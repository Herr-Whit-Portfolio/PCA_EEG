function cSig = EEGPCA(EEG, segm)
%This function is designed to interface with the MATLAB app
%"SingleEEGPCA1_3" - See the same file fore the struct definition of the
%EEG datatype. See also 
cSig = EEG;
%We have to options: 
% 1. The EEG signal may be segmented into equal sequential parts. This
% option would be relevant for examining patterns recorded over an extended
% period.
% 2. The EEG signal may be segmented depending on specific markers set for
% example by experimental software or other instrumentation. Choose this
% if event related potentials are of interest.
if strcmp(segm.type, "sequential")
%     We split the signal into equal segments for each channel.
    cSig.nbseg = double(floor(((EEG.pnts - segm.initOffset) / segm.width)));
    cSig.data = zeros(EEG.nbchan, cSig.nbseg, segm.width);
    cSig.channelName = {EEG.chanlocs.labels};
    for channel = 1:EEG.nbchan
        readStart = segm.initOffset;
        
        for j = 1:(cSig.nbseg - 1)
            cSig.data(channel,j,1:segm.width) = EEG.data(channel, readStart:(readStart + segm.width -1));
            readStart = segm.width + readStart;
        end
        [cSig.weights(channel,:,:), cSig.score(channel,:,:)] = pca(squeeze(cSig.data(channel, :, :)));
        cSig.weights(channel,:,:) = cSig.weights(channel,:,:);
    end
    return
elseif strcmp(segm.type, "eventRelated")
%     If the signal is to be split based on markers, these will be
%     identified in the "label" element and then segmented accordingly.
    cSig.channelName = {EEG.chanlocs.labels};
    cSig.nbseg = 0;
    cSig.swidth = segm.lBound + segm.rBound;
    for i = 1:length(EEG.event)
        for j = 1:length(segm.marker)
            cSig.nbseg = cSig.nbseg + strcmp(EEG.event(i).type, segm.marker{j});
        end
    end
    width = segm.lBound + segm.rBound;
    cSig.data = zeros(EEG.nbchan, cSig.nbseg, width);    
    for channel = 1:EEG.nbchan
        segment = 0;
        for j = 1:length(EEG.event)
            for k = 1:length(segm.marker)
                if max(strcmp(segm.marker{k}, EEG.event(j).type))
                    if(EEG.event(j).latency - segm.lBound) > 0
                        segment = segment + 1;
                        cSig.data(channel, segment, 1:width) = EEG.data(channel,(EEG.event(j).latency - segm.lBound):(EEG.event(j).latency + segm.rBound - 1));
                        
                    end
                end
            end
        end
        [cSig.weights(channel,:,:), cSig.score(channel,:,:)] = pca(squeeze(cSig.data(channel, :, :)));
    end
    return
end
end

function segm = createSegment()
segm.type = "eventRelated";
segm.marker = ["c1", "c2", "c3"];
segm.lBound = 100;
segm.rBound = 200;
end