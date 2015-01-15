function varargout = PrimerDesignTool(varargin)
% PRIMERDESIGNTOOL MATLAB code for PrimerDesignTool.fig
%      PRIMERDESIGNTOOL, by itself, creates a new PRIMERDESIGNTOOL or raises the existing
%      singleton*.
%
%      H = PRIMERDESIGNTOOL returns the handle to a new PRIMERDESIGNTOOL or the handle to
%      the existing singleton*.
%
%      PRIMERDESIGNTOOL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PRIMERDESIGNTOOL.M with the given input arguments.
%
%      PRIMERDESIGNTOOL('Property','Value',...) creates a new PRIMERDESIGNTOOL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PrimerDesignTool_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PrimerDesignTool_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PrimerDesignTool
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PrimerDesignTool_OpeningFcn, ...
                   'gui_OutputFcn',  @PrimerDesignTool_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before PrimerDesignTool is made visible.
function PrimerDesignTool_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PrimerDesignTool (see VARARGIN)

% Choose default command line output for PrimerDesignTool
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PrimerDesignTool wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PrimerDesignTool_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pb1.
function pb1_Callback(hObject, eventdata, handles)
% hObject    handle to pb1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tablestr = get(handles.uitable1,'data');
sizestr = size(tablestr,1);

set(handles.text46,'string', ['No errors']);
set(handles.pb1,'ForegroundColor',[.5 .5 .5]);
set(handles.pb1,'string', 'Evaluating...');
drawnow

newpos = get(handles.text46,'Position');
LinkLength = 0;
for i = 1:sizestr
    if (isempty(tablestr{i,1}) || isempty(tablestr{i,2}))
        sizestr = i-1;
        break
    end
end

set(handles.edit1,'string', sizestr);

for i = 1:sizestr
    seqsPre{i,1} = tablestr{i,1};
    seqsPre{i,2} = tablestr{i,2};
end

Pseq = cell(1,1);
Pseq2 = cell(1,1);

% Association matrix based on http://www.ncbi.nlm.nih.gov/pubmed/19969549
scr_mat = [2 -4 -4 -4; -4 5 -4 -4; -4 -4 5 -4; -4 -4 -4 2]/62*16;

% Number of nucleotides from 3' end to check for association
CATemp = get(handles.edit24,'string');
CAss = str2num(CATemp);

% Initial Minimum Association
tun = 3;

% Minimum primer length
minM = 19;

% Minimum melting temperature
TsetTemp = get(handles.edit25,'string');
Tset = str2num(TsetTemp);

% How much the system weighs association.  Higher is more stringent, but narrows population making larger densities.
Sense = 2; 

% Number of new sequences solved per round.
stseq = 2;

seqs= cell(1,1);

% Homodimerization Parameter
HDimerSt = get(handles.edit23,'string');
HDimer = str2num(HDimerSt);

% Loading user defined sequences in manipulated cell

if size(seqsPre,1) < 2
    stseq = 1;
end

for i = 1:stseq
   
    for j = 1:2
        
        seqs{i,j} = seqsPre{i,j};
        
        if j == 2

            seqs{i,j}=seqreverse(seqcomplement(seqs{i,j}));

        end
        
    end
    
end

% Number of loop/program iterations
SeqMax = (size(seqsPre,1)-stseq)/stseq;

if stseq > 1.5;
    
    if mod(length(seqsPre),2) == 0
        
        odf = 0;
    else
        
        odf = .5;
    
    end
    
else
    
    odf = 0;

end

%///////////////////////////////////////////////////////////////////////
% This is the major loop that calculates best primers, updates primers with best ones, then adds new
% primers to the primer pool
%///////////////////////////////////////////////////////////////////////

for SeqSeq = 1:SeqMax+1+odf
    
    % Adds solve sequences back into main cell, takes reverse
    % complements of reverse sequences then adds new sequences to cell to
    % be solved for.
    if SeqSeq == SeqMax+1+odf && odf == .5
        odf2 = -1;
    else
        odf2 = 0;
    end
    
    if SeqSeq > 1
        
        seqs = Pseq2;
                
        for i = stseq+2*(SeqSeq-1)-1:stseq+2*(SeqSeq-1)+odf2
            
            for j = 1:2

                seqs{i,j}=seqsPre{i,j};
                
                if j == 2
            
                    seqs{i,j}=seqreverse(seqcomplement(seqs{i,j}));
            
                end
    
            end
            
        end
        
    end

    set(handles.text45,'string', num2str(stseq+2*(SeqSeq-1)+odf2));
    drawnow
    pfilt = cell(1,1);
    I = size(seqs);
    Sites = I(1);

    for SC = 1:Sites
        
    % This loop will take a sequence and find primers that melt at 60, do
    % not form hairpins, 40 < GC < 60 and do not have long repeats for all 
    % sequences.  These results will be saved in the pfilt cell.

    % Scheck : Constructed sequence for each loop
    % ForwardP/ReverseP : Appropriate primer sets
    % Pfilt : Cell Storing Primer Sets at row SC
    % SC : Site to be amplified
        reps = 3;
        findex = zeros(1,5);
        rindex = zeros(1,5);
        % Forward Primer Filtering

        seqf = seqs{SC,1};
        N = length(seqf);
        
        for i = 1:N-minM+1

            k = 0;

            for j = 0:N-minM-i+1    

                if k == 0;
                    index = i:i+minM-1+j;
                    Scheck = seqf(index);
                    V = oligoprop(Scheck,'PRIMERCONC',1e-6,'HPLOOP',4);
                    Tm = mean(V.Tm);

                    if Tm > Tset-.5 && Tm < Tset+.5

                        findex(i,1) = length(index);
                        findex(i,2) = V.GC;
                        findex(i,3) = min(size(V.Dimers));
                        findex(i,4) = min(size(V.Hairpins));
                        findex(i,5) = ~cellfun('isempty',regexpi(cellstr(Scheck),'a{3,}|c{3,}|g{3,}|t{3,}','ONCE'));
                        findex(i,6) = ~cellfun('isempty',regexpi(cellstr(Scheck),'a{4,}|c{4,}|g{3,}|t{4,}','ONCE'));
                        findex(i,7) = ~cellfun('isempty',regexpi(cellstr(Scheck),'a{5,}|c{5,}|g{4,}|t{5,}','ONCE'));
                        findex(i,8) = swalign(Scheck,seqcomplement(seqreverse(Scheck)),'SCORINGMATRIX',scr_mat,'GAPOPEN',10,'ALPHA','NT');
                        k = 1;

                    end
                   

                end

            end

        end
        
        if size(findex,2)<7
            B = {['No acceptable primers found for site ',num2str(SC)];[' forward.']}
            B2 = char(B);
            set(handles.text46,'string', B2);
            set(handles.pb1,'string', 'Calculate Primers');
            set(handles.pb1,'ForegroundColor',[0 0 0]);
            drawnow
            error(['No acceptable primers found for site ',num2str(SC),' forward.'])
        end
        
        mix = 0;
        repe = 0;

        % Filters primers and checks if there are appropriate sets.
        while repe == 0

            F1 = find(findex(:,2)<60+mix & findex(:,2)>40-mix & findex(:,8)<HDimer & findex(:,1)<40); % GC Filter
            F2 = find(findex(F1,4)<1 & findex(F1,5)<1); % Hairpin & Repeat Filter
            
            if isempty(F2)
                
                F2 = find(findex(F1,4)<1 & findex(F1,6)<1); % Hairpin & Repeat Filter
                
                if isempty(F2)
                    
                    F2 = find(findex(F1,4)<1 & findex(F1,7)<1); % Hairpin & Repeat Filter            
                    
                end

            end
            
            if isempty(F2)

                mix = mix+2.5;
                if mix > 20
                    set(handles.text46,'string', ['No acceptable primers found for site ',num2str(SC),' forward.']);
                    set(handles.pb1,'string', 'Calculate Primers');
                    set(handles.pb1,'ForegroundColor',[0 0 0]);
                    drawnow
                    error(['No acceptable primers found for site ',num2str(SC),' forward.'])
                end
                
            else
                
                repe = 1;
                
            end
        
        end

        ForwardP = zeros(length(F1(F2)),4);
        ForwardP(:,1) = 1:length(F1(F2)); % Visual Index
        ForwardP(:,2) = F1(F2);           % Start Position (from furthest end)
        ForwardP(:,3) = findex(F1(F2),1); % Length
        ForwardP(:,4) = N-F1(F2)-findex(F1(F2),1)+1; % Length from gen info
 
        pfilt{SC,1} = ForwardP;
        
        % Reverse Primer Filtering

        seqr = seqs{SC,2};
        N = length(seqr); % length of the target sequence
        
        for i = 1:N-minM+1

            k = 0;

            for j = 0:N-minM-i+1    

                if k == 0;

                    index = i:i+minM-1+j;
                    Scheck = seqr(index);
                    V = oligoprop(Scheck,'PRIMERCONC',1e-6,'HPLOOP',4);
                    Tm = mean(V.Tm);

                    if Tm > Tset-.5 && Tm < Tset+.5

                        rindex(i,1) = length(index);
                        rindex(i,2) = V.GC;
                        rindex(i,3) = min(size(V.Dimers));
                        rindex(i,4) = min(size(V.Hairpins));
                        rindex(i,5) = ~cellfun('isempty',regexpi(cellstr(Scheck),'a{3,}|c{3,}|g{3,}|t{3,}','ONCE'));
                        rindex(i,6) = ~cellfun('isempty',regexpi(cellstr(Scheck),'a{4,}|c{4,}|g{3,}|t{4,}','ONCE'));
                        rindex(i,7) = ~cellfun('isempty',regexpi(cellstr(Scheck),'a{5,}|c{5,}|g{4,}|t{5,}','ONCE'));
                        rindex(i,8) = swalign(Scheck,seqcomplement(seqreverse(Scheck)),'SCORINGMATRIX',scr_mat,'GAPOPEN',10,'ALPHA','NT');
                        k = 1;

                    end

                end

            end

        end
        
        if size(rindex,2)<7
            set(handles.text46,'string', ['No acceptable primers found for site ',num2str(SC),' reverse.']);
            set(handles.pb1,'ForegroundColor',[0 0 0]);
            drawnow
            error(['No acceptable primers found for site ',num2str(SC),' reverse.'])
        end
        
        mix = 0;
        repe = 0;
        while repe == 0

            R1 = find(rindex(:,2)<60+mix & rindex(:,2)>40-mix & rindex(:,8)<HDimer & rindex(:,1)<40); % GC Filter
            R2 = find(rindex(R1,4)<1 & rindex(R1,5)<1); % Hairpin & Repeat Filter
            
            if isempty(R2)
                
                R2 = find(rindex(R1,4)<1 & rindex(R1,6)<1); % Hairpin & Repeat Filter
                
                if isempty(R2)
                    
                    R2 = find(rindex(R1,4)<1 & rindex(R1,7)<1); % Hairpin & Repeat Filter            
                    
                end

            end
            
            if isempty(R2)

                mix = mix+2.5;
                if mix > 20
                        set(handles.text46,'string', ['No acceptable primers found for site ',num2str(SC),' reverse.']);
                        set(handles.pb1,'ForegroundColor',[0 0 0]);
                        drawnow
                        error(['No acceptable primers found for site ',num2str(SC),' reverse.'])
                end
                
            else
                
                repe = 1;
                
            end
        
        end

        ReverseP = zeros(length(R1(R2)),4);
        ReverseP(:,1) = 1:length(R1(R2));
        ReverseP(:,2) = R1(R2);
        ReverseP(:,3) = rindex(R1(R2),1);
        ReverseP(:,4) = N-R1(R2)-rindex(R1(R2),1)-1;

        pfilt{SC,2} = ReverseP;

    end

    Sizemat = zeros(SC,2);

    for i = 1:2
    % Calculate number of sites for each primer.  Return warning if the number
    % of sites is small.
        for j = 1:size(pfilt,1)

            Sizemat(j,i) = size(pfilt{j,i},1);

        end

    end
    
    % Beginning association calculations
    % ------------------------------------------------------------
    % ------------------------------------------------------------
    
    % AssocCell contains all interactions with columns corresponding to
    % site + 1 and rows corresponding to site
    SA = 2*Sites-1;
    AssocCell = cell(SA,SA); 
    
    % PFCell contains the possible interactions in 2x2 matrices.
    PFCell = cell(SA,SA);

    % Create a matrix to load information into the cell. rs is a vector
    % which contains the column order which values are evaluated.
    K = repmat(1:SA,SA,1);
    rs = reshape(triu(K)',SA^2,1)';
    rs(rs == 0) = [];

    % Sets up PFCell.
    for i = 1:SA

        count = 0;
        count2 = 0;
        oe = 0;

        for j = 1:SA
            
            % Switches between column 1 or column 2. Count is the vertical
            % index
            if oe == 1;

                PFCell{j,i}(1,:) = [count+1,2];
                oe = 0;

            else

                PFCell{j,i}(1,:) = [count+1,1];
                oe = 1;

            end

            count2 = count2+1;
            
            % Moves index down every two
            if count2 > 1.1

                count = count + 1;
                count2 = 0;
            end

        end

    end


    for i = 1:SA

        count = 0;
        count2 = 1;
        oe = 0;

        for j = 1:SA
            if oe == 1;

                PFCell{i,j}(2,:) = [count+1,1];
                oe = 0;

            else

                PFCell{i,j}(2,:) = [count+1,2];
                oe = 1;

            end

            count2 = count2+1;
            if count2 > 1.1

                count = count + 1;
                count2 = 0;
            end

        end

    end

    imax = SA; % Number of sites to evaluate before shifting rows.
    j = 1; % Current Row
    
    % for loop is 1 to Number of interacting primer types (sites we have 6 interactions)
    
    for i = 1:(2*Sites^2-Sites)
    % This loop will calculate association matricies for all possible
    % primer pairs.  For N sites we have 2N^2-N sites.
    % 
    % These lines determine where the association matrix is loaded.

        in = rs(i);
        
        % Looking up primer indicies to evaluate.
        Set1 = PFCell{j,in}(1,:);
        Set2 = PFCell{j,in}(2,:);

        % Finding primer positions, length to evaluate.
        Filt1 = pfilt{Set1(1),Set1(2)};
        Filt2 = pfilt{Set2(1),Set2(2)};
        
        % Sequences to evaluate
        Seq1 = seqs{Set1(1),Set1(2)};
        
        Seq2 = seqs{Set2(1),Set2(2)};
        
        % Number of different primers
        Sf1 = size(Filt1,1);
        Sf2 = size(Filt2,1);
        
        % Temporary Association Matrix
        Assoc = zeros(Sf1(1),Sf2(1));
        
        % Reconstruct sequences from length data and sequence data.
        for w = 1:Sf1
            % First set of sequences
            index = Filt1(w,2):(Filt1(w,2)+Filt1(w,3)-1);
            Fs = Seq1(index);

            for x = 1:Sf2
                
                % Second set of sequences, and calculates association.
                index = Filt2(x,2):Filt2(x,2)+Filt2(x,3)-1;
                Rs = Seq2(index);
                W3(1) = swalign(Fs(length(Fs)-CAss:length(Fs)),seqcomplement(seqreverse(Rs)),'SCORINGMATRIX',scr_mat,'EXTENDGAP',10,'GAPOPEN',10,'ALPHA','NT');
                Assoc(w,x) = W3(1);

            end

        end   

        % Temp data to permenant stoarge
        AssocCell{j,in} = Assoc;
        
        % Switches to a new row
        if i == imax
            imax = imax+SA-j;
            j = j+1;
        end

    end

    % Now have AssocCell, which is a cell with lookup values for all primer
    % associations.  pfilt contains all the sequence information
    % corresponding to the indicies of the matrix is AssocCell.
    
    % Threshold and loop controllers.
    LoopControl2 = 0;
    intsort = 0;
    
    while LoopControl2 == 0
    % This loop will calculate least interacting sites.  This is the
    % difficult part of the program because combinatorial statistics causes
    % this problem to become huge. If no sites can be found at a certain
    % condition then threshold (tun) is increased by 1 and the problem
    % repeats itself.
        Set1 = 0;

        % Initial set satisfying association conditions.
        [CA RA] = find(AssocCell{1,1}<=1.1+tun);
        
        % If no associations match criteria, threshold is increased.  If
        % there are interactions Set1 is generated.  Set1 contains all
        % possible combinations and is a conserved variable through the
        % loop.  Set2 feeds into Set1.  Rows for Set 1 are different
        % primers, and columns are different sets
        
        if isempty(CA)
            LoopControl1 = 1;
        else
            Set1 = zeros(2*Sites,length(CA));
            Set1(1,1:length(CA)) = CA';
            Set1(2,:) = RA';
            LoopControl1 = 0;
        end
        
        for k = 2:SA
        % k is the column in the Association matrix we loop at.  Loop
        % starts at k = 2, because k = 1 is previously calculated.
            intsort = intsort + 1;
            if LoopControl1 == 0
            % If values exists for Set1 under current threshold conditions.
                kmax = 0;
                for m = 1:k
                % m is the row in the Association matrix. This is the new
                % set of primers that we test for each loop.
                    if k > kmax
                        kmax = k;
                    end
                    
                    if LoopControl1 == 0
                        
                        % Find appropriate associations at m,k
                        [CA RA] = find(AssocCell{m,k} <= 1.1+tun);
                        
                        if isempty(CA)
                            % No associations found, need larger
                            % threshold
                            LoopControl1 = 1;

                        else
                            
                            % Potential Primers to keep, add a column on
                            % Set1 then calculate acceptable sites to keep
                            CAU2 = unique(CA);
                            settrk = cell(1,1);
                            sizeset = 0;
                            for L = 1:length(CAU2)
                            
                                settrk{L} = find(Set1(m,:) == CAU2(L));
                                sizeset = size(settrk{L},2)*sum(CAU2(L) == CA)+sizeset;
                                
                            end
                            
                            Set2 = zeros(2*Sites,sizeset);
                            LC = 0; %Load Counter
                            for L = 1:size(settrk,2) % Index in settrk which is an array containing the indicies of indentified CAU2 sites
                                
                                for P = 1:length(settrk{L}) % P is CAU2 sites of L in Set1
                                    
                                    for Q = 1:sum(CAU2(L) == CA) % Q is the number of incidents of a CAU2 site with respect to association
                                                            
                                        LC = LC+1;
                                        Set2(:,LC) = Set1(:,settrk{L}(P));
                                        TempR = RA(CAU2(L) == CA);
                                        Set2(k+1,LC)= TempR(Q);
                                        
                                    end
                                    
                                end
                                
                            end
                            
                            % if no sites found then increase threshold
                            if isempty(Set2)
                                LoopControl1 = 1;
                                % If no site for this primer agrees
                                % with the pool then increase
                                % threshold.

                            else

                                Set1 = unique(Set2','rows')';

                            end

                        end

                    end

                end
                
            end
            
        end
        
        if LoopControl1 == 0
            
            LoopControl2 = 1;
            
        else
            
            tun = tun+.5;
            
        end
       
        
    end

    % This section ranks the results of Set 1 based on similarities.
    for w = 1:30

        S1 = size(Set1,2);
        %Matrix which store association values.
        ACheck = 0;

        for r = 1:S1

            ind = 0;

            for i = 1:SA

                for j = i:SA

                    ind = ind+1;
                    ACheck(ind,r) = AssocCell{i,j}(Set1(i,r),Set1(j+1,r));

                end

            end

        end
        
        I = find(sum(ACheck == 16-w/2,1) == min(sum(ACheck == 16-w/2,1)));
        Set1 = Set1(:,I);

    end

    % This section ranks the results of Set 1 based length from information

    SP = size(pfilt,1);
    ACheck = 0;
    ACheck2 = 0;
    SumL = 0;

    for i = 1:size(Set1,2)

        Ind = 0;

        for j = 1:SP

            for k = 1:2
                Ind = Ind+1;
                AA = pfilt{j,k}(Set1(Ind,i),:);
                ACheck(i,Ind) = AA(4);
                ACheck2(i,Ind) = AA(3);

            end

        end

    end

    [I J] = min(sum(ACheck,2));
    Best = J(1);
    SumL = sum(ACheck(Best,:),2);
    LinkLength = LinkLength+I(1)+SumL;

    % List Primers

    Ind = 0;
    seqs2 = cell(1,1);
    PSeq = cell(1,1);

    for j = 1:SP

        for k = 1:2
            
            seqs2{j,k} = seqs{j,k};
            Ind = Ind+1;
            AA = pfilt{j,k}(Set1(Ind,Best),:);
            Pseq{j,k} = seqs2{j,k}(AA(2):AA(2)+AA(3)-1);
        end

    end
    Pseq2 = Pseq;
    set(handles.uitable2,'data', Pseq);
    tun3 = tun+1;
    AData2 = num2str(tun3);
    set(handles.text44,'string', AData2);
    drawnow
    
end

tun3 = tun+1;

AData2 = num2str(tun3);

set(handles.text44,'string', AData2);
set(handles.pb1,'ForegroundColor',[0 0 0]);
set(handles.pb1,'string', 'Calculate Primers');

function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.pushbutton2,'ForegroundColor',[.5 .5 .5]);
set(handles.pushbutton2,'string', 'Evaluating...');
SCd = get(handles.text45,'string');
tuntemp = get(handles.edit8,'string');
SC = str2num(SCd);
disp(SCd)
SP = SC;
scr_mat = [2 -4 -4 -4; -4 5 -4 -4; -4 -4 5 -4; -4 -4 -4 2]/62*16;
Off = 0;
set(handles.text47,'string', num2str(Off));
findex = 0;
Trk = 0;
LLeft = 0:SC-1;
RandCounter = 0;
tun2 = str2num(tuntemp);
TsetTemp = get(handles.edit27,'string');
Tset = str2num(TsetTemp)-10;
Pseq = get(handles.uitable2,'data');
Pseq2 = Pseq;
drawnow
SeqLTemp = get(handles.edit26,'string');
SeqL = str2num(SeqLTemp);

tl = 1.5;

% Homodimerization Parameter
HDimerSt = get(handles.edit23,'string');
HDimer = str2num(HDimerSt);
% Max tuning interaction parameter between population and strand of
% interest

for i = 1:size(Pseq,1)
    for j = 1:2
    maxtun(i,j) = tun2;
        Fseq = Pseq{i,j};
        for i2 = 1:size(Pseq,1)
            for j2 = 1:2
                    Rseq = Pseq{i2,j2};
                    [I,J] = swalign(Fseq,seqcomplement(seqreverse(Rseq)),'SCORINGMATRIX',scr_mat,'GAPOPEN',10,'ALPHA','NT');
                    if I+1 > maxtun(i,j)
                        maxtun(i,j) = I+1;
                    end
            end
        end
    end
end

maxtun

maxtun2 = [maxtun(2:end,1), maxtun(1:end-1,2)];
mmax = max(maxtun2')
umax = unique(mmax)
% Solving for the least interacting primers.
maxtracker = 1;
Icount = 1;

while Off < SC-1 && RandCounter < 1e10
    
    I = find(umax(maxtracker) == mmax);
    lI = length(I);
    
    if Icount > lI
        maxtracker = maxtracker + 1;
        Icount = 1;
        I = find(umax(maxtracker) == mmax)
    end
    
    RandCounter = RandCounter + 1;
    Hi1 = cell(1,1);
    HS1 = 0;
    
    % Makes a random linker sequence
    AssocPrimer = 0;
    Rs = randseq(SeqL);
    V = oligoprop(Rs,'PRIMERCONC',1e-6,'HPLOOP',4);
    Mt = mean(V.Tm);

    if abs(Tset+10-Mt)<1  % Checks to see that the linker sequence is the temperature we want
       findex(1,1) = V.GC;
       findex(1,2) = ~cellfun('isempty',{V.Hairpins}');
       findex(1,3) = ~cellfun('isempty',regexpi(cellstr(Rs),'a{4,}|c{4,}|g{4,}|t{4,}','ONCE'));
       findex(1,4) = swalign(Rs,seqreverse(seqcomplement(Rs)),'SCORINGMATRIX',scr_mat,'GAPOPEN',10,'ALPHA','NT');

       if findex(1,1)<60 && findex(1,1)>40 && findex(1,2)<1 && findex(1,3)<1 && findex(1,4)<tun2-tl% GC, Hairpin, Repeat Filter, Self-Association Dimer
           %Performing association analysis with primer set.
           
            set(handles.text48,'string', num2str(RandCounter));
            drawnow
            for j = 1:SP

                for k = 1:2
                    STemp = Pseq2{j,k};
                    LTemp = length(STemp);
                    W3 = swalign(STemp,seqcomplement(seqreverse(Rs)),'SCORINGMATRIX',scr_mat,'GAPOPEN',10,'ALPHA','NT');
                    W4 = swalign(STemp(LTemp-12:LTemp),seqcomplement(seqreverse(Rs)),'SCORINGMATRIX',scr_mat,'GAPOPEN',10,'ALPHA','NT');
                    if W3 > tun2-tl && W4 > tun2-tl-1
                        AssocPrimer = 1;
                    end

                end

            end

            if AssocPrimer == 0

                Rs2 = seqreverse(seqcomplement(Rs));

                for j = 1:SP

                    for k = 1:2
                        STemp = Pseq2{j,k};
                        LTemp = length(STemp);
                        W3 = swalign(STemp,seqcomplement(seqreverse(Rs2)),'SCORINGMATRIX',scr_mat,'GAPOPEN',10,'ALPHA','NT');
                        W4 = swalign(STemp(LTemp-12:LTemp),seqcomplement(seqreverse(Rs2)),'SCORINGMATRIX',scr_mat,'GAPOPEN',10,'ALPHA','NT');
                        if W3 > tun2-tl && W4 > tun2-tl-1
                               AssocPrimer = 1;
                        end

                    end

                end

            end

            if AssocPrimer == 0
                HS1 = 0;
                HS2 = 0;
                HS3 = 0;
                HS4 = 0;
                HS5 = 0;
                HS6 = 0;
                ToCheck = setdiff(LLeft,Trk);
                TCL = 0;

                for j = ToCheck
                    TCL = TCL+1;
                    TS1 = [Rs2,Pseq{j,2}];
                    TS2 = [Rs,Pseq{j+1,1}];
                    VTemp1 = oligoprop(TS1,'PRIMERCONC',1e-6,'HPLOOP',4);
                    VTemp2 = oligoprop(TS2,'PRIMERCONC',1e-6,'HPLOOP',4);
                    HS1(TCL) = (size(VTemp1.Hairpins,1) < .1);
                    HS2(TCL) = (size(VTemp2.Hairpins,1) < .1);
                    HS3(TCL) = swalign(TS1,seqcomplement(seqreverse(TS1)),'SCORINGMATRIX',scr_mat,'GAPOPEN',10,'ALPHA','NT')<=HDimer;
                    HS4(TCL) = swalign(TS2,seqcomplement(seqreverse(TS2)),'SCORINGMATRIX',scr_mat,'GAPOPEN',10,'ALPHA','NT')<=HDimer;
                    HS5(TCL) = (~cellfun('isempty',regexpi(cellstr(TS1),'a{5,}|c{5,}|g{4,}|t{5,}','ONCE'))<1);
                    HS6(TCL) = (~cellfun('isempty',regexpi(cellstr(TS2),'a{5,}|c{5,}|g{4,}|t{5,}','ONCE'))<1);
                end

                HPv = find(HS1 & HS2 & HS3 & HS4 & HS5 & HS6 == 1);
                HPv = ToCheck(HPv);
        
                if isempty(HPv) || not(isempty(intersect(HPv,Trk)))
                    AssocPrimer = 1;
                else
                    
                    %Checking the association between possible new primers
                    %and the population
                   HPvfc = 0;
                   HPvrc = 0;
                   HPvf = 0;
                   HPvr = 0;
                                   
                   for i = 1:length(HPv)
                        W3 = 2;
                        for j = 1:SP

                            for k = 1:2

                                TS1 = [Rs2,Pseq{HPv(i),2}];
                                W3temp = swalign(Pseq2{j,k},seqcomplement(seqreverse(TS1)),'SCORINGMATRIX',scr_mat,'GAPOPEN',10,'ALPHA','NT');
                                if W3temp > W3
                                    W3 = W3temp;
                                end

                            end

                        end
                        
                        if W3 <= maxtun(HPv(i),2)
                            HPvfc = HPvfc + 1;
                            HPvf(HPvfc) = HPv(i);
                        end

                    end


                   for i = 1:length(HPv)

                        W3 = 2;
                        for j = 1:SP

                            for k = 1:2
                                TS2 = [Rs,Pseq{HPv(i)+1,1}];
                                W3temp = swalign(Pseq2{j,k},seqcomplement(seqreverse(TS2)),'SCORINGMATRIX',scr_mat,'GAPOPEN',10,'ALPHA','NT');
                                if W3temp > W3
                                    W3 = W3temp;
                                end

                            end

                        end
                        
                        if W3 <= maxtun(HPv(i)+1,1)

                               HPvrc = HPvrc + 1;
                               HPvr(HPvrc) = HPv(i);

                        end

                   end

                    HPv = intersect(HPvf,HPvr);
                    set(handles.text48,'string', num2str(RandCounter));
                    drawnow
                    if isempty(HPv) || sum(HPv) == 0
                        AssocPrimer = 1;
                    end

                    if AssocPrimer == 0;

                        HI = setdiff(HPv,Trk);
                        HI2 = intersect(HI,I);
                        
                        if not(isempty(HI2))
                            
                            Trk(length(Trk)+1) = HI2(1);
                            Off = Off+1;
                            Pseq2{HI2(1),2} = [Rs2,Pseq{HI2(1),2}];
                            Pseq2{HI2(1)+1,1} = [Rs,Pseq{HI2(1)+1,1}];
                            Icount = Icount + 1;
                            set(handles.text47,'string', num2str(Off));
                            set(handles.uitable3,'data',Pseq2);
                            drawnow
                        
                        end

                    end


                end

            end

       end

    end

end

set(handles.text48,'string', num2str(RandCounter));
set(handles.uitable3,'data', Pseq2);

set(handles.pushbutton7,'ForegroundColor',[0 0 0]);
set(handles.pushbutton2,'ForegroundColor',[0 0 0]);
set(handles.pushbutton2,'string', 'Calculate Linkers');


function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on edit3 and none of its controls.
function edit3_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over edit5.
function edit5_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FN,PN] = uigetfile('*.csv','Select the csv');
FN2 = [PN FN];
fid = fopen(FN2);
KK = textscan(fid,'%s%s%s','Delimiter',',');
KN = [KK{2},KK{3}];
set(handles.uitable1,'data', KN);
SK = size(KK{1},1);
SK2 = 10-SK;

for i = 1:SK2
     KK{1}(10-i+1) = cellstr(' ');
end

set(handles.uitable1, 'rowname', KK{1});
set(handles.uitable2, 'rowname', KK{1});
set(handles.uitable3, 'rowname', KK{1});

function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Si = get(handles.edit1,'string');
Si = str2num(Si);
M1 = get(handles.uitable3, 'rowname');
M1 = M1(1:Si);
M2 = get(handles.uitable3,'data');
Si2 = get(handles.text47,'string');
Si2 = str2num(Si2);
if size(M2,1) == Si2+1
    
    M = [M1 M2];
    [I J] = uiputfile('*.csv','Save primers and linkers');
    path = [J I];

    fid = fopen(path, 'w');
    
    for row=1:Si
        fprintf(fid, '%s,%s,%s\n', M{row,:});
    end
    
    fclose
    
else
   
    
    
end
    


% --- Executes when figure1 is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit23_Callback(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit23 as text
%        str2double(get(hObject,'String')) returns contents of edit23 as a double


% --- Executes during object creation, after setting all properties.
function edit23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit24_Callback(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit24 as text
%        str2double(get(hObject,'String')) returns contents of edit24 as a double


% --- Executes during object creation, after setting all properties.
function edit24_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit25_Callback(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit25 as text
%        str2double(get(hObject,'String')) returns contents of edit25 as a double


% --- Executes during object creation, after setting all properties.
function edit25_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit26_Callback(hObject, eventdata, handles)
% hObject    handle to edit26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit26 as text
%        str2double(get(hObject,'String')) returns contents of edit26 as a double


% --- Executes during object creation, after setting all properties.
function edit26_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit27_Callback(hObject, eventdata, handles)
% hObject    handle to edit27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit27 as text
%        str2double(get(hObject,'String')) returns contents of edit27 as a double


% --- Executes during object creation, after setting all properties.
function edit27_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Si = get(handles.edit1,'string');
Si = str2num(Si);
M1 = get(handles.uitable2, 'rowname');
M1 = M1(1:Si);
M2 = get(handles.uitable2,'data');
Si2 = get(handles.text47,'string');
Si2 = str2num(Si2);

    M = [M1 M2];
    [I J] = uiputfile('*.csv','Save primers and linkers');
    path = [J I];

    fid = fopen(path, 'w');
    
    for row=1:Si
        fprintf(fid, '%s,%s,%s\n', M{row,:});
    end
    
    fclose
