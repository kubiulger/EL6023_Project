function [antsym,antind] = hPRGPrecodeTEST(siz,nstartgrid,portsym,portind,F)
    
    % Get the number of precoder resource groups NPRG from the precoder
    % array F
    NPRG = size(F,3);
    
    % Get the number of resource blocks from the size vector
    NRB = siz(1) / 12;
    
    % Get the PRG numbers (1-based) for each CRB in the whole carrier
    prgset = getPRGSet(NRB,nstartgrid,NPRG);
    
    % Establish the dimensionality of the grid of port indices
    portsiz = siz;
    ndims = max(length(portsiz),3);
    portsiz(end) = size(F,1);
    
    % Calculate 1-based RE subscripts from port indices
    [subs{1:ndims}] = ind2sub(portsiz,portind);
    resubs = subs{1};
    
    % Calculate 0-based CRB subscripts from RE subscripts
    crbsubs = floor((resubs - 1) / 12);
    
    % Calculate 1-based PRG subscripts from CRB subscripts
    prgsubs = prgset(crbsubs + 1);

    % Perform precoding to produce antenna symbols and antenna indices
    [antsym,antind] = hPrecode(siz,portsym,portind,F,prgsubs);

end

% Get PRG numbers (1-based) for a carrier with NRB resource blocks,
% starting CRB 'nstartgrid' and NPRG precoder resource groups
function prgset = getPRGSet(NRB,nstartgrid,NPRG)

    Pd_BWP = ceil((NRB + nstartgrid) / NPRG);
    prgset = repmat(1:NPRG,[Pd_BWP 1]);
    prgset = reshape(prgset(nstartgrid + (1:NRB).'),[],1);

end

% Precoding of symbols and projection of the corresponding indices
function [symout,indout] = hPrecode(siz,symin,indin,F,prgsubs)
    
    % dimensionality information
    ndims = max(length(siz),3);
    symindims = arrayfun(@(x)size(symin,x),1:(ndims-1));
    inplanedims = symindims(2:end);
    outplanedim = size(F,2);
    NPRG = size(F,3);
    outdims = [siz(1:(ndims-1)) outplanedim];
    
    % calculate 1-based port subscripts and unique set of ports from
    % port indices
    [subs{1:ndims}] = ind2sub(siz,indin);
    lastinplanesubs = subs{end};
    lastinplanes = unique(lastinplanesubs(:)).';
    
    % create output antenna indices and empty output antenna symbols
    [symout,indout] = nrExtractResources(indin,zeros(outdims,'like',symin));
    
    % for each port
    for plane = lastinplanes
        
        % port indices / symbols which correspond to this port
        thisplane = (lastinplanesubs==plane);
        
        % for each PRG
        for prg = 1:NPRG

            % port indices / symbols which correspond to this PRG
            thisprg = (prgsubs==prg);

            % beamform the symbols for this port and PRG using the
            % appropriate beam, producing the symbols for all antennas (for
            % this port and PRG), note that logical indexing of 'symin'
            % loses its shape (restored below if required)
            bfsym = symin(thisplane & thisprg) * F(plane,:,prg);

            % if there are beamformed symbols for this PRG
            if (~isempty(bfsym))
            
                % if port symbols / indices have multiple planes, restore
                % the correct shape for the beamformed symbols - this
                % assumes that the symbols / indices are in multiple
                % dimensions with symbols / indices sorted according to
                % increasing linear indices
                if (numel(inplanedims)>1)
                    bfdims = {[] inplanedims(1:end-1) outplanedim};
                    bfsym = reshape(bfsym,bfdims{:});
                end

                % calculate rows of antenna indices that correspond to the
                % current PRG
                outplaneprgrows = all(reshape(thisprg,size(thisprg,1),[]),2);

                % add the beamformed symbols to the antenna symbols for
                % this PRG (may overlap with beamformed symbols for other
                % ports in the same PRG)
                if (numel(inplanedims)>1)
                    outplaneprgsubs = {find(outplaneprgrows)};
                    outplaneprgsubs(2:(1+numel(inplanedims))) = {':'};
                    symout(outplaneprgsubs{:}) = symout(outplaneprgsubs{:}) + bfsym;
                else
                    symout(outplaneprgrows,:) = symout(outplaneprgrows,:) + bfsym;
                end
                
            end
            
        end
        
    end
    
end


