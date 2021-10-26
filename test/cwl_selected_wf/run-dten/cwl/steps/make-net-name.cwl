label: make-net-name
id: make-net-name
cwlVersion: v1.0
class: ExpressionTool

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
outputs:
  net-name:
    type: string

expression: |
    ${
      var ns = inputs.netpre.replace(/\s/g, '')
      var net = ns+"_beta"+inputs.beta +"_mu"+inputs.mu+"_w"+inputs.w;
      return { "net-name": net };
     }
