## Project 2: Speech Enhancement

file structure:


    │  evaculate_stoi.m calculate STOI for given mode
    │  main.m this is the initial file, can see the structure of the whole project
    │  main_all_noise.m get the result of all noise
    │  main_noise1.m get the result only with noise 1
    │  main_noise2.m get the result only with noise 2
    │  main_noise3.m get the result only with noise 3
    │  main_noise4.m get the result only with noise 4
    │  plot_STOI.m plot the STOI result for different beamformer
    │
    ├─data
    │      aritificial_nonstat_noise.wav
    │      babble_noise.wav
    │      clean_speech.wav
    │      clean_speech_2.wav
    │      impulse_responses.mat
    │      Speech_shaped_noise.wav
    │
    ├─data_result result saved in .mat file
    │
    ├─imgs images saved in .png file
    │
    └─util
            GEVD.m carry out GEVD algorithm
            speechnormalize.m normalize speech signal to [-1,1]
            stoi.m