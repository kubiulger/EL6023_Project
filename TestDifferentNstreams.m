load roomPathData.mat

LinkStates = pathData.linkState;
LOS = find(LinkStates==1); %indices for LOS links
NLOS = find(LinkStates==2); %Used later

rxNum = LOS(randi([1,length(LOS)])); %Selecting an RX number from LOS indices
txPos = pathData.txPos; %Getting tx Position
rxPos = pathData.rxPos(rxNum,:); %Getting rx Position

npaths = pathData.npaths(rxNum,:); %Number of paths
gain = pathData.gain(rxNum,:); %Gain of each path
dly = pathData.dly(rxNum,:); %Delay of each path

% Angles 
aoaAz = pathData.aoaAz(rxNum,:)';
aoaEl = pathData.aoaEl(rxNum,:)';
aodAz = pathData.aodAz(rxNum,:)';
aodEl = pathData.aodEl(rxNum,:)';

%We may wish to Rx and Tx towards each other so calculating some parameters
trVec = rxPos' - txPos; %Vector from tx to rx
trVec = trVec/norm(trVec);
zvec = [0,0,1]'; %Antenna elements are directed towards +z by default

tiltAx = cross(zvec,[trVec;0]);
tiltAngle = 90; %Since tiltAx,trVec and zVec are mutually orthogonal

fc = pathData.fc; %Carrier Frequency
elemTx = design(patchMicrostrip, fc); %Simple antenna element
elemTx.Tilt = 90; 
elemTx.TiltAxis = tiltAx';
elemRx = design(patchMicrostrip, fc); %Simple antenna element
elemRx.Tilt = 90; 
elemRx.TiltAxis = -tiltAx';

nantgNB = [4,4]; % 16 Antennas at transmitter
nantUE = 8; % 8 Antennas at receiver

lambda = physconst('LightSpeed')/fc; %Wavelength
dstep = lambda/2; %Distance between antennas

arrgNB = phased.URA(nantgNB,dstep,"ArrayNormal","x"); %URA for transmitter
arrUE = phased.ULA(nantUE,dstep,"ArrayAxis","y"); %ULA for receiver

arrPlatformgNB = ArrayPlatform('elem', elemTx, 'arr', arrgNB, 'fc', fc);
arrPlatformgNB.computeNormMatrix();

arrPlatformUE =  ArrayPlatform('elem', elemRx, 'arr', arrUE, 'fc', fc);
arrPlatformUE.computeNormMatrix();

SubcarrierSpacing = 120;    % SCS in kHZ
NRB = 61;  % number of resource blocks
nscPerRB = 12;  % number of sub-carriers per RB

carrierConfig = nrCarrierConfig(...
    'NSizeGrid', NRB, 'SubcarrierSpacing', SubcarrierSpacing);
waveformConfig = nrOFDMInfo(carrierConfig);

fdchan = FDMIMOChan(carrierConfig, 'txArrPlatform', arrPlatformgNB, 'rxArrPlatform', arrPlatformUE, ...
    'aoaAz', aoaAz, 'aodAz', aodAz, 'aoaEl', aoaEl, 'aodEl', aodEl,  ...
    'gain', gain, 'dly', dly, 'fc', fc);

frameNum = 0;
slotNum = 0;
[chanGrid, Q] = fdchan.step(frameNum, slotNum);

mcsInd = [10,11,12,13,14,15,16];

Modulation = '16QAM';
target = [340,378,434,790,553,616,658]/1024;
err = zeros(length(target),length(nStreamsVec));

for j=1:length(target)
    targetCodeRate = target(j);
    pdschConfig = nrPDSCHConfig(...
        'Modulation', Modulation, ...
        'PRBSet', (0:NRB-1), ...
        'SymbolAllocation', [1, waveformConfig.SymbolsPerSlot-1], ...
        'PTRS', nrPDSCHPTRSConfig(),...
        'NumLayers', 1);


    nStreamsVec = [1,2,3,4,5,6]; %Chosen for dmeonstration

    for i=1:length(nStreamsVec)
        errSubToT=0;
        totSubTot=0;
        for k=1:10
            tx = NRgNBTxFD(carrierConfig, pdschConfig, 'targetCodeRate', targetCodeRate);
            nStreams=nStreamsVec(i);
            [txGrid,Fprecoder] = tx.step(nStreams,Q);
            G = pagemtimes(chanGrid,Fprecoder); %Pre-coded channel, it's assumed that RX knows this.
            AvgSNR = 15; %Let's assume this is the SNR we have, we will pick noise variance to match this.
            noiseVar = sum(mean(abs(G).^2,[3,4]),'all')/db2pow(AvgSNR)/nantUE/nStreams;

            %Now we pass txGrid through the channel and add noise to obtain received
            %grid txGrid
            rxGridNoNoise = pagemtimes(chanGrid,txGrid);
            rxGrid = rxGridNoNoise + sqrt(noiseVar/2)*(randn(size(rxGridNoNoise))+1j*randn(size(rxGridNoNoise)));
            rx = NRUERxFD(carrierConfig, pdschConfig, 'targetCodeRate', targetCodeRate);
            rx.step(rxGrid, G, noiseVar,nStreams); %Demodulation

            rxBits = rx.rxBits;
            txBits = tx.txBits;
            ErrorBits = length(find(rxBits ~= txBits));
            TotalBits = length(txBits);
            errSubTot = errSubToT + ErrorBits;
            totSubTot = totSubTot + TotalBits;
        end
        err(j,i)=errSubTot/totSubTot;       
    end
end


Sim = struct;
Sim.rxNum = rxNum;
Sim.Error = err;
Sim.MCS = mcsInd;
Sim.NumLayers = nStreamsVec;
save('Simulation','Sim')
