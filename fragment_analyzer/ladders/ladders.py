import numpy as np

LADDERS = {
    "LIZ": {
        "sizes": np.array([
            80,
            100,
            114,
            120,
            140,
            160,
            180,
            200,
            214,
            220,
            240,
            260,
            280,
            300,
            314,
            320,
            340,
            360,
            380,
            400,
            414,
            420,
            440,
            460,
            480,
            500,
            514,
            520,
            540,
            560,
            580,
            600
        ]),
        "distance": 30, 
        "height": 100, 
        "max_ladder_trace_distance": 300
    },
    "ROX": {
        "sizes": np.array([
            79,
            90,
            105,
            131,
            151,
            182,
            201,
            254,
            306,
            337,
            362,
            425,
            486,
            509,
            560,
            598,
            674,
            739,
            799,
            902,
            1007,
        ]),
        "distance": 30, 
        "height": 100, 
        "max_ladder_trace_distance": 300
    }
}


CHANNELS = {"LIZ": "DATA205", "ROX": "DATA12"}
