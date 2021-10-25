{
    "class": "CommandLineTool",
    "requirements": [
        {
            "class": "InlineJavascriptRequirement"
        }
    ],
    "inputs": [
        {
            "type": "Any",
            "default": {
                "baz": "zab1",
                "b az": 2,
                "b'az": true,
                "b\"az": null,
                "buz": [
                    "a",
                    "b",
                    "c"
                ]
            },
            "id": "#main/bar"
        }
    ],
    "outputs": [
        {
            "id": "#params_inc.yml/t1",
            "type": "Any",
            "outputBinding": {
                "outputEval": "$(inputs)"
            }
        },
        {
            "id": "#params_inc.yml/t2",
            "type": "Any",
            "outputBinding": {
                "outputEval": "$(inputs.bar)"
            }
        },
        {
            "id": "#params_inc.yml/t3",
            "type": "Any",
            "outputBinding": {
                "outputEval": "$(inputs['bar'])"
            }
        },
        {
            "id": "#params_inc.yml/t4",
            "type": "Any",
            "outputBinding": {
                "outputEval": "$(inputs[\"bar\"])"
            }
        },
        {
            "id": "#params_inc.yml/t5",
            "type": "Any",
            "outputBinding": {
                "outputEval": "$(inputs.bar.baz)"
            }
        },
        {
            "id": "#params_inc.yml/t6",
            "type": "Any",
            "outputBinding": {
                "outputEval": "$(inputs['bar'].baz)"
            }
        },
        {
            "id": "#params_inc.yml/t7",
            "type": "Any",
            "outputBinding": {
                "outputEval": "$(inputs['bar'][\"baz\"])"
            }
        },
        {
            "id": "#params_inc.yml/t8",
            "type": "Any",
            "outputBinding": {
                "outputEval": "$(inputs.bar['baz'])"
            }
        },
        {
            "id": "#params_inc.yml/t9",
            "type": "Any",
            "outputBinding": {
                "outputEval": "$(inputs.bar['b az'])"
            }
        },
        {
            "id": "#params_inc.yml/t10",
            "type": "Any",
            "outputBinding": {
                "outputEval": "$(inputs.bar['b\\'az'])"
            }
        },
        {
            "id": "#params_inc.yml/t11",
            "type": "Any",
            "outputBinding": {
                "outputEval": "$(inputs.bar[\"b'az\"])"
            }
        },
        {
            "id": "#params_inc.yml/t12",
            "type": "null",
            "outputBinding": {
                "outputEval": "$(inputs.bar['b\"az'])"
            }
        },
        {
            "id": "#params_inc.yml/t13",
            "type": "Any",
            "outputBinding": {
                "outputEval": "-$(inputs.bar.baz)"
            }
        },
        {
            "id": "#params_inc.yml/t14",
            "type": "Any",
            "outputBinding": {
                "outputEval": "-$(inputs['bar'].baz)"
            }
        },
        {
            "id": "#params_inc.yml/t15",
            "type": "Any",
            "outputBinding": {
                "outputEval": "-$(inputs['bar'][\"baz\"])"
            }
        },
        {
            "id": "#params_inc.yml/t16",
            "type": "Any",
            "outputBinding": {
                "outputEval": "-$(inputs.bar['baz'])"
            }
        },
        {
            "id": "#params_inc.yml/t17",
            "type": "Any",
            "outputBinding": {
                "outputEval": "$(inputs.bar.baz) $(inputs.bar.baz)"
            }
        },
        {
            "id": "#params_inc.yml/t18",
            "type": "Any",
            "outputBinding": {
                "outputEval": "$(inputs['bar'].baz) $(inputs['bar'].baz)"
            }
        },
        {
            "id": "#params_inc.yml/t19",
            "type": "Any",
            "outputBinding": {
                "outputEval": "$(inputs['bar'][\"baz\"]) $(inputs['bar'][\"baz\"])"
            }
        },
        {
            "id": "#params_inc.yml/t20",
            "type": "Any",
            "outputBinding": {
                "outputEval": "$(inputs.bar['baz']) $(inputs.bar['baz'])"
            }
        },
        {
            "id": "#params_inc.yml/t21",
            "type": "Any",
            "outputBinding": {
                "outputEval": "$(inputs.bar['b az']) $(inputs.bar['b az'])"
            }
        },
        {
            "id": "#params_inc.yml/t22",
            "type": "Any",
            "outputBinding": {
                "outputEval": "$(inputs.bar['b\\'az']) $(inputs.bar['b\\'az'])"
            }
        },
        {
            "id": "#params_inc.yml/t23",
            "type": "Any",
            "outputBinding": {
                "outputEval": "$(inputs.bar[\"b'az\"]) $(inputs.bar[\"b'az\"])"
            }
        },
        {
            "id": "#params_inc.yml/t24",
            "type": "Any",
            "outputBinding": {
                "outputEval": "$(inputs.bar['b\"az']) $(inputs.bar['b\"az'])"
            }
        },
        {
            "id": "#params_inc.yml/t25",
            "type": "Any",
            "outputBinding": {
                "outputEval": "$(inputs.bar.buz[1])"
            }
        },
        {
            "id": "#params_inc.yml/t26",
            "type": "Any",
            "outputBinding": {
                "outputEval": "$(inputs.bar.buz[1]) $(inputs.bar.buz[1])"
            }
        },
        {
            "id": "#params_inc.yml/t27",
            "type": "null",
            "outputBinding": {
                "outputEval": "$(null)"
            }
        },
        {
            "id": "#params_inc.yml/t28",
            "type": "int",
            "outputBinding": {
                "outputEval": "$(inputs.bar.buz.length)"
            }
        }
    ],
    "baseCommand": "true",
    "id": "#main",
    "cwlVersion": "v1.2"
}