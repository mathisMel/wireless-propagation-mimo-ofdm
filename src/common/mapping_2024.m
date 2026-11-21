function [symb_tx] = mapping_2024(bit_tx,Nbps,modulation, nb_source)

% INPUTS:
% - bit_tx : vector of input bits 
% - Nbps : number of bits per symbol
% - modulation : 'pam' or 'qam'
%
% OUTPUTS:
% - symb_tx : vector of ouput symbols (variance 1)

Nsymb = size(bit_tx,1)/Nbps; % Number of symbols

bit_tx2 = permute(reshape(bit_tx,Nbps,Nsymb,nb_source), [2 1 3]);
int_tx = zeros(Nsymb,nb_source);
int_txI = zeros(Nsymb,nb_source);
int_txQ = zeros(Nsymb,nb_source);
switch modulation,
    
    case 'bpsk'
        
        % Gray to binary
        mapp_tx(:,1,:) = bit_tx2(:,1,:);
        for ii = 2:Nbps,
           mapp_tx(:,ii,:) = xor( mapp_tx(:,ii-1,:) , bit_tx2(:,ii,:) ); 
        end

        % Binary to integer
        for i = 1:nb_source
            int_tx(:,i) = bi2de(fliplr(mapp_tx(:,:,i)));
        end

        % Integer to symbol
        sigma = sqrt(sum(([0:2^Nbps-1]-(2^Nbps-1)/2).^2)/2^Nbps); 
        symb_tx = 1/sigma * (int_tx - (2^Nbps-1)/2);

    case 'qam'
        
        % REAL PART
        NbpsI = Nbps/2;
        bit_tx2I = bit_tx2(:,1:NbpsI,:);
        
        % Gray to binary
        mapp_txI(:,1,:) = bit_tx2I(:,1,:);
        for ii = 2:NbpsI,
           mapp_txI(:,ii,:) = xor( mapp_txI(:,ii-1,:) , bit_tx2I(:,ii,:) ); 
        end

        % Binary to integer
        for i = 1:nb_source
            int_txI(:,i) = bi2de(fliplr(mapp_txI(:,:,i)));
        end
        % int_txI(:,1) = bi2de(fliplr(mapp_txI(:,:,1)));
        % int_txI(:,2) = bi2de(fliplr(mapp_txI(:,:,2)));

        % Integer to symbol
        sigmaI = sqrt(sum(([0:2^NbpsI-1]-(2^NbpsI-1)/2).^2)/2^NbpsI); 
        symb_txI = 1/sigmaI/sqrt(2) * (int_txI - (2^NbpsI-1)/2);
        
        
        % IMAGINARY PART
        NbpsQ = Nbps/2;
        bit_tx2Q = bit_tx2(:,NbpsQ+1:end,:);
        
        % Gray to binary
        mapp_txQ(:,1,:) = bit_tx2Q(:,1,:);
        for ii = 2:NbpsQ,
           mapp_txQ(:,ii,:) = xor( mapp_txQ(:,ii-1,:) , bit_tx2Q(:,ii,:) ); 
        end

        % Binary to integer
        for i = 1:nb_source
            int_txQ(:,i) = bi2de(fliplr(mapp_txQ(:,:,i)));
        end
        % int_txQ(:,1) = bi2de(fliplr(mapp_txQ(:,:,1)));
        % int_txQ(:,2) = bi2de(fliplr(mapp_txQ(:,:,2)));

        % Integer to symbol
        sigmaQ = sqrt(sum(([0:2^NbpsQ-1]-(2^NbpsQ-1)/2).^2)/2^NbpsQ); 
        symb_txQ = 1/sigmaQ/sqrt(2) * (int_txQ - (2^NbpsQ-1)/2);
       
        
        % COMPLEX SYMBOL
        symb_tx = symb_txI + j*symb_txQ;
       
end