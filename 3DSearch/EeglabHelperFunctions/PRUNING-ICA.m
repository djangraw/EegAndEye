% Subject 2 ICA components judged to be bad (at a quick glance)
% Targets
EEG = pop_subcomp( EEG, [1   2   3   4   5   6  25  32  43  44  49  50  51  52  59  61  62  64], 0);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 4,'setname','3DS-2-targsac-ica-pruned','overwrite','on','gui','off'); 
%Distractors
EEG = pop_subcomp( EEG, [1   2   3   4   5   6  15  27  32  33  36  37  39  41  43  46  55  56  60  61  63], 0);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 5,'setname','3DS-2-distsac-ica-pruned','overwrite','on','gui','off'); 