function [OFDMBlocks, dataIndices,pilotIndices,zeroIndices] = CreateOFDMBlocks_ELEC422_423(mappedData,mappedPilots,pilotPercentage,totalSC,num0SC,numOFDMBlocks, nb_source)

% inserting the pilots and data into OFDM subcarriers
% INPUTS:
%  - mappedData: mapped symbols of the data bits
%  - mappedPilots: mapped symbols of the pilot bits
%  - pilotPercentage: percentage of nonzero subcarriers that are pilots
%  - totalSC: total number of SC
%  - num0SC: number of SC that are zero
%  - numOFDMBlocks: amount of OFDM blocks 
% OUTPUTS: 
%  - OFDMBlocks: all OFDM blocks, rows are the different blocks and columns
%                the diff SC
%  - dataIndices: the indices of data SC
%  - pilotIndices: the indices of pilot SC
%  - zeroIndices: the indices of zero SC

    numNon0SC = 1/2*(totalSC - 2*num0SC - 1 - 3); % 3 center frequency 0's, one more on the right
    
    dataBlockSize = ceil(1/pilotPercentage)-1;
    pilotBlockSize = ceil(1/(1-pilotPercentage))-1;
    
    numFullBlocks = floor(2*numNon0SC/(dataBlockSize+pilotBlockSize));
    numAtEnd = 2*numNon0SC - numFullBlocks*(dataBlockSize+pilotBlockSize);

    pilotIndices = (1:2*numNon0SC).';
    %mappedPilots = mappedPilots(1:numOFDMBlocks);
    
    if numAtEnd == dataBlockSize
        %in this case, the last data block is fully there and not the pilots
        totalDataSymbols = (numFullBlocks+1)*dataBlockSize;
        totalPilotSymbols = 2*numNon0SC - totalDataSymbols;
    
        numTotalData = totalDataSymbols*numOFDMBlocks;
        dataSymbols = mappedData(1:numTotalData,:);
    
        %find indices of pilots, WERKT
        pilotIndices(end-dataBlockSize+1:end,:) = []; %remove last data block indices
        pilotIndices = reshape(pilotIndices,dataBlockSize+pilotBlockSize,[]);
        pilotIndices(1:dataBlockSize,:) = []; %remove the data indices
   
    elseif numAtEnd == 0
        %in this case, each row has integer amount of blocks of data and pilots
        totalDataSymbols = (numFullBlocks)*dataBlockSize;
        totalPilotSymbols = 2*numNon0SC - totalDataSymbols;
    
        numTotalData = totalDataSymbols*numOFDMBlocks;
        dataSymbols = mappedData(1:numTotalData,:);
    
        %find indices of pilots, WERKT
        pilotIndices = reshape(pilotIndices,dataBlockSize+pilotBlockSize,[]);
        pilotIndices(1:dataBlockSize,:) = []; %remove the data indices

    elseif numAtEnd < dataBlockSize
        %in this case, the data is partially there
        totalDataSymbols = (numFullBlocks)*dataBlockSize + numAtEnd;
        totalPilotSymbols = 2*numNon0SC - totalDataSymbols;
    
        numTotalData = totalDataSymbols*numOFDMBlocks;
        dataSymbols = mappedData(1:numTotalData,:);
    
        %find indices of pilots, WERKT
        pilotIndices(end-numAtEnd+1:end,:) = []; %remove partial data block indices
        pilotIndices = reshape(pilotIndices,dataBlockSize+pilotBlockSize,[]);
        pilotIndices(1:dataBlockSize,:) = [];
    
    else
        %in this case, the pilots are partially there
        totalDataSymbols = (numFullBlocks+1)*dataBlockSize;
        totalPilotSymbols = 2*numNon0SC - totalDataSymbols;
    
        numTotalData = totalDataSymbols*numOFDMBlocks;
        dataSymbols = mappedData(1:numTotalData,:);

        %for the pilots: remove the pilot indices for the ones added to the end
        correctPilots = pilotIndices(end-numAtEnd+dataBlockSize+1:end); %grab the pilots blocks that are partially there
        pilotIndices(end-numAtEnd+1:end,:) = []; %remove these indices
        
        %find indices of pilots, WERKT
        pilotIndices = reshape(pilotIndices,pilotBlockSize+dataBlockSize,[]);
        pilotIndices(1:dataBlockSize,:) = [];
        pilotIndices = reshape(pilotIndices,[],1);
        pilotIndices = [pilotIndices; correctPilots];
        sort(pilotIndices);

    end
    
    dataIndices = setdiff((1:2*numNon0SC).',pilotIndices);
    %dataSymbols = permute(dataSymbols, [2 1]);
    dataSymbols = permute(reshape(dataSymbols,totalDataSymbols,numOFDMBlocks,nb_source), [2 1 3]);
    mappedPilots = permute(repmat(mappedPilots, 1, 1, totalPilotSymbols),[1, 3, 2]);


    totalVector = zeros(numOFDMBlocks ,2*numNon0SC,nb_source);
    totalVector(:,dataIndices,:) = dataSymbols;
    totalVector(:,pilotIndices,:) = mappedPilots;
    OFDMBlocks = totalVector; %output OFDM blocks
%     preamble = repmat(sign(randn(1, 2*numNon0SC)), 2,1);
%     OFDMBlocks(1:2,:) = preamble;
    [dataIndices,pilotIndices,zeroIndices] = ShiftIndices(dataIndices,pilotIndices,num0SC,totalSC);
    
end



function [dataIndices,pilotIndices,zeroIndices] = ShiftIndices(dataIndices,pilotIndices,num0SC,totalSC)
    % shift the indices with the zero carriers
    % INPUTS:
    %  - data- & pilotIndices
    %  - number of 0 SC and total number of SC
    % OUTPUTS:
    %  - shifted data- & pilotIndices 
    %  - zeroIndices

    numNon0SC = 1/2*(totalSC - 2*num0SC - 1 - 3); %3 center frequency 0's, one more on the right

    DCIndices = numNon0SC+1:numNon0SC+3; %remove DC subcarriers from indexes
    negIndices = find(dataIndices<=numNon0SC); %shift half of pilot and data indexes right with 3
    negIndices = reshape(negIndices,1,[]);
    posIndices = find(dataIndices>numNon0SC);
    posIndices = reshape(posIndices,1,[]);
    dataIndices = [reshape(dataIndices(negIndices),1,[]), reshape(dataIndices(posIndices)+3,1,[])];
    dataIndices = setdiff(dataIndices,DCIndices);
    
    negIndices = find(pilotIndices<=numNon0SC);
    negIndices = reshape(negIndices,1,[]);
    posIndices = find(pilotIndices>numNon0SC);
    posIndices = reshape(posIndices,1,[]);
    pilotIndices = [reshape(pilotIndices(negIndices),1,[]), reshape(pilotIndices(posIndices)+3,1,[])];
    pilotIndices = setdiff(pilotIndices,DCIndices);
    
    dataIndices(dataIndices>2*numNon0SC+3) = [];%remove indices too high
    pilotIndices(pilotIndices>2*numNon0SC+3) = [];
    
    zeroIndices = [1:num0SC,numNon0SC+num0SC+1:numNon0SC+num0SC+3,totalSC-num0SC:totalSC];

end