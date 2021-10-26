#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: ExpressionTool
id: make-net-name
label: make-net-name
requirements:
  InlineJavascriptRequirement: {}
inputs:
  beta:
    type: double
  mu:
    type: double
  w:
    type: double
  netpre:
    type: string
expression: |
  ${
    var ns = inputs.netpre.replace(/\s/g, '')
    var net = ns+"_beta"+inputs.beta +"_mu"+inputs.mu+"_w"+inputs.w;
    return { "net-name": net };
   }
outputs:
  net-name:
    type: string
