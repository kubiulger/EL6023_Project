classdef NRUERxFD < matlab.System
    % 5G NR UR receiver class implemented in frequency domain
    properties
        % Configuration
        carrierConfig;   % Carrier configuration
        pdschConfig;     % Default PDSCH config
        waveformConfig;  % Waveform config
        
        % OFDM grid
 %       rxGrid;
        
        % Transport block data for last transmission
        targetCodeRate = 490/1024;  % Target code rate
        trBlkSizes;                 % Transport block size
        
        % Received data in last slots
        pdschEq;       % Equalized PDSCH symbols
        rxBits;        % RX bits
        
        % DLSCH decoder
        decDLSCH;
        
        
    end
    methods
        function obj = NRUERxFD(carrierConfig, pdschConfig, ...
                varargin)
            % Constructor
            
            % Save the carrier and PDSCH configuration
            obj.carrierConfig = carrierConfig;
            obj.pdschConfig = pdschConfig;
            
            % Create the waveform configuration from the carrier
            % configuration
            obj.waveformConfig = nrOFDMInfo(obj.carrierConfig);
            
            % Set parameters from constructor arguments
            if nargin >= 1
                obj.set(varargin{:});
            end
            
            % Create DLSCH decoder
            obj.decDLSCH = nrDLSCHDecoder('MultipleHARQProcesses', false, ...
                'TargetCodeRate', obj.targetCodeRate, ...
                'LDPCDecodingAlgorithm', 'Layered belief propagation');
            
        end
    end
    methods (Access = protected)
        
        
        function stepImpl(obj, rxGridInput, G, noiseVar, nStreams)
            %Hermitian of the precoded channel (known at the receiver)
            Gconj = permute(conj(G),[2,1,3,4]);
            Flmmse = zeros(size(Gconj));
            Qlmmse = zeros([nStreams,nStreams,732,14]);
            alp = 1/noiseVar/16;
            for i=1:732
                for j = 1:14
                    %LMMSE decoder at each resource element
                    Qlmmse(:,:,i,j) =(alp*(Gconj(:,:,i,j)*G(:,:,i,j))...
                        +eye(nStreams))^-1;
                    Flmmse(:,:,i,j)=alp*Qlmmse(:,:,i,j)*Gconj(:,:,i,j); 
                end
            end
            %Mean of Q is taken so that we get an idea on average SNR
            Qlmmse = mean(Qlmmse,[3,4]);
            %So now we calculate equalized 1/SNR which is going to be a
            %real number but we still take the real to make sure
            Qii = real(diag(Qlmmse));
            eqNoise = (1./(1./Qii-1));
            
            %Equalized channel, it should be close to identity matrix at high SNR
            Heq = pagemtimes(Flmmse,G);
            %Performing Equalizer on the received grid
            rxGrid3D = pagemtimes(Flmmse,rxGridInput);
           

            for i=1:nStreams
                % Demodulates and decodes one slot of data
                chanGrid = squeeze(Heq(i,i,:,:));
                rxGrid = squeeze(rxGrid3D(i,:,:));
                noiseVar1 = eqNoise(i);
                % Get PDSCH received symbols and channel estimates
                % from received grid 
                [pdschInd,pdschInfo] = nrPDSCHIndices(...
                    obj.carrierConfig, obj.pdschConfig);
                [pdschRx, pdschHest] = nrExtractResources(pdschInd, rxGrid,...
                    chanGrid);

                % TODO:  Perform the MMSE equalization using the
                % nrEqualizeMMSE() function.
                % Use the PDSCH Rx symbols, PDSCH channel estimate and noise
                % variance as the input.  Store the equalized symbols in
                % obj.pdschEq and  channel state information in a structure,
                % csi.            
                %    [obj.pdschEq,csi] = nrEqualizeMMSE(...);
                [obj.pdschEq,csi] = nrEqualizeMMSE(pdschRx,pdschHest,noiseVar1);

                % TODO:  Get the LLRs with the nrPDSCHDecode() function.
                % Use carrier and PDSCH configuration, the equalized symbols,
                % and the noise variance, noiseVar.
                %    [obj.pdschEq,csi] = nrEqualizeMMSE(...);   
                [dlschLLRs,rxSym] = nrPDSCHDecode(obj.carrierConfig,obj.pdschConfig,obj.pdschEq,noiseVar);
                % Scale LLRs by EbN0.  
                % The csi value computed in the nrEqualizeMMSE()
                % function is csi = |pdschHest|^2 + noiseVar.
                % Also, the Eb/N0 = snrEq/Qm where Qm is the number of bits 
                % per symbol and snrEq is the SNR after equalization,
                %
                %   snrEq = (|pdschHest|^2 + noiseVar)/noiseVar = csi/noiseVar
                %
                % Hence, Eb/N0 = csi/(noiseVar*Qm).
                % Since the LLRs from the nrPDSCHDecode function are 
                % already scaled by 1/noiseVar, we multiply them by  csi/Qm.
                csi = nrLayerDemap(csi); % CSI layer demapping
                numCW = length(csi);
                for cwIdx = 1:numCW
                    Qm = length(dlschLLRs{cwIdx})/length(rxSym{cwIdx}); % bits per symbol
                    csi{cwIdx} = repmat(csi{cwIdx}.',Qm,1);   % expand by each bit per symbol
                    dlschLLRs{cwIdx} = dlschLLRs{cwIdx} .* csi{cwIdx}(:);   % scale
                end

                % Compute the extra overhead from the PT-RS
                Xoh_PDSCH = 6*obj.pdschConfig.EnablePTRS;

                % Calculate the transport block size based on the PDSCH
                % allocation and target code rate
                obj.trBlkSizes = nrTBS(obj.pdschConfig.Modulation,obj.pdschConfig.NumLayers,...
                    numel(obj.pdschConfig.PRBSet),pdschInfo.NREPerPRB,...
                    obj.targetCodeRate,Xoh_PDSCH);
                obj.decDLSCH.TransportBlockLength = obj.trBlkSizes;

                % Reset the soft buffer
                harqId = 0;
                obj.decDLSCH.resetSoftBuffer(harqId);

                % TODO:  Decode the bits with the obj.decDLSCH() method.
                % Use the scaled LLRs from above. Use a redundancy version, 
                % rv = 0,  since we are not using HARQ in this lab.
                %    rv = 0;
                %    obj.rxBit = obj.decDLSCH(...);
    %             rv = [0,0];
                rv=0;
                %Get bits in this layer
                newBits = obj.decDLSCH(dlschLLRs,obj.pdschConfig.Modulation,obj.pdschConfig.NumLayers,rv);
                %Add it to rest
                obj.rxBits =[obj.rxBits; newBits];
            end
        end
        
    end
end

