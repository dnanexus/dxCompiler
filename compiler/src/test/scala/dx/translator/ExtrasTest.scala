package dx.translator

import com.dnanexus.exceptions.ResourceNotFoundException
import dx.Assumptions.isLoggedIn
import dx.Tags.{ApiTest, EdgeTest}
import dx.PermissionDeniedException
import dx.api.{DxAccessLevel, DxApi, DxPath, DxProject}
import dx.core.ir.Value
import dx.core.ir.Value.VString
import dx.translator.ExtrasJsonProtocol._
import dx.util.Logger
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import spray.json._

class ExtrasTest extends AnyFlatSpec with Matchers {
  private val dxApi = DxApi()(Logger.Quiet)
  private val testProject: String = "dxCompiler_playground"
  private lazy val project: DxProject = {
    try {
      dxApi.resolveProject(testProject)
    } catch {
      case _: ResourceNotFoundException =>
        throw new Exception(
            s"""|Could not find project ${testProject}, you probably need to be logged into
                |the platform on staging.""".stripMargin
        )
    }
  }

  private def getIdFromName(name: String): String = {
    dxApi.resolveDataObject(name, Some(project)).id
  }

  private lazy val appletId: String = getIdFromName("/release_test/mummer_nucmer_aligner")
  private lazy val fileId: String = dxApi.resolveDataObject("Readme.md", Some(project)).id

  it should "recognize restartable entry points" in {
    val runtimeAttrs: JsValue =
      """|{
         |  "defaultTaskDxAttributes" : {
         |     "runSpec" : {
         |       "restartableEntryPoints": "all"
         |     }
         |  }
         |}""".stripMargin.parseJson

    val extras = Extras.parse(runtimeAttrs)
    extras.defaultTaskDxAttributes should be(
        Some(DxAppJson(Some(DxRunSpec(None, None, Some("all"), None))))
    )
  }

  it should "invalid runSpec I" in {
    val ex1 =
      """|{
         | "defaultTaskDxAttributes": {
         |     "timeoutPolicy": {
         |        "main": {
         |          "hours": 12
         |        }
         |     }
         |  }
         |}
         |""".stripMargin

    assertThrows[DeserializationException] {
      Extras.parse(ex1.toJson)
    }
  }

  it should "invalid runSpec II" in {
    val ex2 =
      """|{
         | "defaultTaskDxAttributes" : {
         |   "runSpec": {
         |     "timeoutPolicy": {
         |        "main": {
         |          "hours": 12
         |        }
         |     }
         |   }
         |  }
         |}
         |""".stripMargin

    assertThrows[Exception] {
      Extras.parse(ex2.parseJson)
    }
  }

  it should "invalid runSpec III" in {
    val ex3 =
      """|{
         | "defaultTaskDxAttributes" : {
         |   "runSpec": {
         |     "access" : {
         |        "project": "CONTRIBUTE__XD",
         |        "allProjects": "INVAL",
         |        "developer": true
         |     }
         |   }
         |  }
         |}
         |""".stripMargin

    assertThrows[Exception] {
      Extras.parse(ex3.parseJson)
    }
  }

  it should "parse the run spec" in {
    val runSpecValid =
      """|{
         | "defaultTaskDxAttributes" : {
         |   "runSpec": {
         |     "executionPolicy": {
         |        "restartOn": {
         |          "*": 3
         |        }
         |     },
         |     "timeoutPolicy": {
         |        "*": {
         |          "hours": 12
         |        }
         |     },
         |     "access" : {
         |        "project": "CONTRIBUTE",
         |        "allProjects": "VIEW",
         |        "network": [
         |           "*"
         |         ],
         |        "developer": true
         |     }
         |   }
         |  }
         |}
         |""".stripMargin

    val js = runSpecValid.parseJson
    val extras = Extras.parse(js)
    extras.defaultTaskDxAttributes should be(
        Some(
            DxAppJson(
                Some(
                    DxRunSpec(
                        Some(
                            DxAccess(Some(Vector("*")),
                                     Some(DxAccessLevel.Contribute),
                                     Some(DxAccessLevel.View),
                                     Some(true),
                                     None)
                        ),
                        Some(DxExecPolicy(Some(Map("*" -> 3)), None)),
                        None,
                        Some(DxTimeout(None, Some(12), None))
                    )
                ),
                None
            )
        )
    )
  }

  it should "parse complex execution policy" in {
    val runSpec =
      """|{
         | "defaultTaskDxAttributes" : {
         |   "runSpec": {
         |     "executionPolicy": {
         |        "restartOn": {
         |           "UnresponsiveWorker": 2,
         |           "JMInternalError": 0,
         |           "ExecutionError": 4
         |        },
         |        "maxRestarts" : 5
         |     }
         |    }
         |  }
         |}
         |""".stripMargin

    val js = runSpec.parseJson
    val extras = Extras.parse(js)

    val restartPolicy: Map[String, Long] =
      Map("UnresponsiveWorker" -> 2, "JMInternalError" -> 0, "ExecutionError" -> 4)
    extras.defaultTaskDxAttributes should be(
        Some(
            DxAppJson(
                Some(DxRunSpec(None, Some(DxExecPolicy(Some(restartPolicy), Some(5))), None, None)),
                None
            )
        )
    )
  }

  it should "recognize error in complex execution policy" in {
    val runSpec =
      """|{
         | "defaultTaskDxAttributes" : {
         |   "runSpec": {
         |     "executionPolicy": {
         |        "restartOn": {
         |           "UnresponsiveWorker_ZZZ": 2,
         |           "ExecutionError": 4
         |        },
         |        "maxRestarts" : 5
         |     }
         |    }
         |  }
         |}
         |""".stripMargin

    val js = runSpec.parseJson
    assertThrows[Exception] {
      Extras.parse(js)
    }
  }

  // Need to include method to upload the README.md file.
  ignore should "parse the custom_reorg object" taggedAs ApiTest in {
    assume(isLoggedIn)

    // inputs is Readme.md file in
    val reorg: JsValue =
      s"""|{
          | "customReorgAttributes" : {
          |    "appUri" :"${appletId}",
          |    "configFile" : "${fileId}"
          |  }
          |}
          |""".stripMargin.parseJson

    val extras = Extras.parse(reorg)
    extras.customReorgAttributes shouldBe Some(CustomReorgSettings(appletId, Some(fileId)))
  }

  it should "throw DeserializationException due to missing applet id" taggedAs ApiTest in {
    assume(isLoggedIn)
    val inputs: String = "dx://file-123456"
    val reorg: JsValue =
      s"""|{
          | "customReorgAttributes" : {
          |    "configFile" : "${inputs}"
          |  }
          |}
          |""".stripMargin.parseJson

    val thrown = intercept[DeserializationException] {
      Extras.parse(reorg)
    }

    thrown.getMessage should be("appUri must be specified in the customReorgAttributes section")
  }

  it should "throw DeserializationException due to missing inputs in custom_reorg section" taggedAs ApiTest in {
    assume(isLoggedIn)

    val reorg: JsValue =
      s"""|{
          | "customReorgAttributes" : {
          |    "appUri" : "${appletId}"
          |  }
          |}
          |""".stripMargin.parseJson

    val thrown = intercept[DeserializationException] {
      Extras.parse(reorg)
    }

    thrown.getMessage should include(
        "'configFile' must be specified as a valid DNAnexus"
    )
  }

  it should "Allow inputs to be null in custom reorg" taggedAs ApiTest in {
    assume(isLoggedIn)
    val reorg: JsValue =
      s"""|{
          | "customReorgAttributes" : {
          |    "appUri" : "${appletId}",
          |    "configFile" : null
          |  }
          |}
          |""".stripMargin.parseJson

    val extras = Extras.parse(reorg)
    extras.customReorgAttributes shouldBe Some(CustomReorgSettings(appletId))
  }

  it should "throw DeserializationException due to invalid applet ID" taggedAs ApiTest in {
    assume(isLoggedIn)
    // invalid applet ID
    val appId: String = "applet-123456"
    val reorg: JsValue =
      s"""|{
          |  "customReorgAttributes" : {
          |      "appUri" : "${appId}",
          |      "configFile": null
          |   }
          |}
          |""".stripMargin.parseJson

    val thrown = intercept[DeserializationException] {
      Extras.parse(reorg)
    }

    thrown.getMessage should be(
        s"Invalid reorg app(let) ID applet-123456"
    )
  }

  ignore should "throw ResourceNotFoundException due to non-existent applet" taggedAs ApiTest in {
    assume(isLoggedIn)
    // non-existent (made up) applet ID
    val appletId: String = "applet-mPX7K2j0Gv2K2jXF75Bf21v2"
    val reorg: JsValue =
      s"""|{
          |  "customReorgAttributes" : {
          |      "appUri" : "${appletId}",
          |      "configFile": null
          |   }
          |}
          |""".stripMargin.parseJson
    val thrown = intercept[Exception] {
      Extras.parse(reorg)
    }

    thrown.getMessage should be(
        s"""Error running command dx describe $appletId --json"""
    )

  }
  it should "throw DeserializationException due to invalid file ID" taggedAs ApiTest in {
    assume(isLoggedIn)
    val inputs: String = "file-1223445"
    val reorg: JsValue =
      s"""|{
          |  "customReorgAttributes" : {
          |      "appUri" : "${appletId}",
          |      "configFile": "${inputs}"
          |   }
          |}
          |""".stripMargin.parseJson

    val thrown = intercept[DeserializationException] {
      Extras.parse(reorg)
    }
    thrown.getMessage should include(
        "'configFile' must be specified as a valid DNAnexus"
    )
  }
  it should "throw ResourceNotFoundException due to non-existant file" taggedAs ApiTest in {
    assume(isLoggedIn)
    val inputs: String = "dx://file-AZBYlBQ0jy1qpqJz17gpXFf8"
    val reorg: JsValue =
      s"""|{
          |  "customReorgAttributes" : {
          |      "appUri" : "${appletId}",
          |      "configFile": "${inputs}"
          |   }
          |}
          |""".stripMargin.parseJson

    val thrown = intercept[ResourceNotFoundException] {
      Extras.parse(reorg)
    }

    val fileId: String = inputs.replace(DxPath.DxUriPrefix, "")

    thrown.getMessage should be(s""""${fileId}" is not a recognized ID""")
  }

  ignore should "throw PermissionDeniedException due to applet not having contribute access in the project" taggedAs ApiTest in {
    assume(isLoggedIn)
    val appletId: String = getIdFromName("/release_test/Sum ")
    val reorg: JsValue =
      s"""|{
          | "customReorgAttributes" : {
          |    "appUri" : "${appletId}",
          |    "configFile": "null"
          |  }
          |}
          |""".stripMargin.parseJson

    val thrown = intercept[PermissionDeniedException] {
      Extras.parse(reorg)
    }

    thrown.getMessage should be(
        s"ERROR: App(let) for custom reorg stage ${appletId} does not " +
          s"have CONTRIBUTOR or ADMINISTRATOR access and this is required."
    )
  }

  it should "take app id as well as applet id for custom reorg" taggedAs (ApiTest, EdgeTest) in {
    assume(isLoggedIn)
    val appId: String = dxApi.resolveApp("cloud_workstation").id
    val reorg: JsValue =
      s"""|{
          |  "customReorgAttributes" : {
          |    "appUri" : "${appId}",
          |    "configFile": null
          |  }
          |}
          |""".stripMargin.parseJson

    val extras = Extras.parse(reorg)
    extras.customReorgAttributes shouldBe Some(CustomReorgSettings(appId))
  }

  it should "generate valid JSON execution policy" in {
    val expectedJs: JsValue =
      """|{
         |  "restartOn": {
         |    "*": 5
         |  },
         |  "maxRestarts" : 4
         |}
         |""".stripMargin.parseJson

    val execPolicy = DxExecPolicy(Some(Map("*" -> 5)), Some(4))
    execPolicy.toJson shouldBe expectedJs
  }

  it should "generate valid JSON timeout policy" in {
    val expectedJs: JsValue =
      """|{
         |  "timeoutPolicy": {
         |    "*": {
         |      "hours": 12,
         |      "minutes" : 30
         |    }
         |  }
         |}
         |""".stripMargin.parseJson
    val runSpec = DxRunSpec(None, None, None, Some(DxTimeout(None, Some(12), Some(30))))
    runSpec.toJson shouldBe expectedJs
  }

  it should "handle default runtime attributes" in {
    val runtimeAttrs: JsValue =
      """|{
         |  "defaultRuntimeAttributes" : {
         |    "docker" : "quay.io/encode-dcc/atac-seq-pipeline:v1"
         |  }
         |}""".stripMargin.parseJson

    val extras = Extras.parse(runtimeAttrs)
    val dockerOpt: Option[Value] = extras.defaultRuntimeAttributes.flatMap(_.get("docker"))
    dockerOpt match {
      case None =>
        throw new Exception("Wrong type for dockerOpt")
      case Some(docker) =>
        docker should equal(VString("quay.io/encode-dcc/atac-seq-pipeline:v1"))
    }
  }

  it should "handle invalid default runtime attributes" in {
    val rtInvalid: JsValue =
      """|{
         | "defaultRuntimeAttributes" : {
         |     "preemptible": "0",
         |     "bootDiskSizeGb": "10"
         | }
         |}""".stripMargin.parseJson

    assertThrows[DeserializationException] {
      Extras.parse(rtInvalid)
    }
  }

  it should "accept per task attributes" in {
    val runSpec: JsValue =
      """|{
         | "defaultTaskDxAttributes" : {
         |   "runSpec": {
         |     "timeoutPolicy": {
         |        "*": {
         |          "hours": 12
         |        }
         |     }
         |   }
         |  },
         | "perTaskDxAttributes" : {
         |   "Add": {
         |      "runSpec": {
         |        "timeoutPolicy": {
         |          "*": {
         |             "minutes": 30
         |          }
         |        }
         |      }
         |    },
         |    "Multiply" : {
         |      "runSpec": {
         |        "timeoutPolicy": {
         |          "*": {
         |            "minutes": 30
         |          }
         |        },
         |        "access" : {
         |          "project": "UPLOAD"
         |        }
         |      }
         |    }
         |  }
         |}
         |""".stripMargin.parseJson

    val extras = Extras.parse(runSpec)
    extras.defaultTaskDxAttributes should be(
        Some(
            DxAppJson(Some(
                          DxRunSpec(
                              None,
                              None,
                              None,
                              Some(DxTimeout(None, Some(12), None))
                          )
                      ),
                      None)
        )
    )
    extras.perTaskDxAttributes should be(
        Some(
            Map(
                "Multiply" -> DxAppJson(
                    Some(
                        DxRunSpec(
                            Some(DxAccess(None, Some(DxAccessLevel.Upload), None, None, None)),
                            None,
                            None,
                            Some(DxTimeout(None, None, Some(30)))
                        )
                    ),
                    None
                ),
                "Add" -> DxAppJson(
                    Some(DxRunSpec(None, None, None, Some(DxTimeout(None, None, Some(30))))),
                    None
                )
            )
        )
    )
  }

  it should "include optional details and runSpec in per task attributes" in {
    val runSpec: JsValue =
      """|{
         | "defaultTaskDxAttributes" : {
         |   "runSpec": {
         |     "timeoutPolicy": {
         |        "*": {
         |          "hours": 12
         |        }
         |     }
         |   }
         |  },
         | "perTaskDxAttributes" : {
         |   "Add": {
         |      "runSpec": {
         |        "timeoutPolicy": {
         |          "*": {
         |             "minutes": 30
         |          }
         |        }
         |      },
         |      "details": {
         |        "upstreamProjects": [
         |          {
         |            "name": "GATK4",
         |            "repoUrl": "https://github.com/broadinstitute/gatk",
         |            "version": "GATK-4.0.1.2",
         |            "license": "BSD-3-Clause",
         |            "licenseUrl": "https://github.com/broadinstitute/LICENSE.TXT",
         |            "author": "Broad Institute"
         |          }
         |        ]
         |      }
         |    },
         |    "Multiply" : {
         |      "runSpec": {
         |        "timeoutPolicy": {
         |          "*": {
         |            "minutes": 30
         |          }
         |        },
         |        "access" : {
         |          "project": "UPLOAD"
         |        }
         |      }
         |    }
         |  }
         |}
         |""".stripMargin.parseJson

    val extras = Extras.parse(runSpec)
    extras.defaultTaskDxAttributes should be(
        Some(
            DxAppJson(Some(
                          DxRunSpec(
                              None,
                              None,
                              None,
                              Some(DxTimeout(None, Some(12), None))
                          )
                      ),
                      None)
        )
    )

    extras.perTaskDxAttributes should be(
        Some(
            Map(
                "Add" -> DxAppJson(
                    Some(DxRunSpec(None, None, None, Some(DxTimeout(None, None, Some(30))))),
                    Some(
                        DxDetails(
                            Some(
                                Vector(
                                    DxLicense("GATK4",
                                              "https://github.com/broadinstitute/gatk",
                                              "GATK-4.0.1.2",
                                              "BSD-3-Clause",
                                              "https://github.com/broadinstitute/LICENSE.TXT",
                                              "Broad Institute")
                                )
                            )
                        )
                    )
                ),
                "Multiply" -> DxAppJson(
                    Some(
                        DxRunSpec(
                            Some(DxAccess(None, Some(DxAccessLevel.Upload), None, None, None)),
                            None,
                            None,
                            Some(DxTimeout(None, None, Some(30)))
                        )
                    ),
                    None
                )
            )
        )
    )
  }

  it should "accept treeTurnaroundTimeThreshold in (default) tasks and workflows" in {
    val runSpec: JsValue =
      """|{
         | "defaultTaskDxAttributes": {
         |   "treeTurnaroundTimeThreshold": 1
         |  },
         | "perTaskDxAttributes": {
         |   "Add": {
         |      "treeTurnaroundTimeThreshold": 2
         |    },
         |    "Multiply" : {
         |      "treeTurnaroundTimeThreshold": 3
         |    }
         |  },
         | "defaultWorkflowDxAttributes": {
         |   "treeTurnaroundTimeThreshold": 5
         | },
         | "perWorkflowDxAttributes": {
         |   "NestedWf": {
         |      "treeTurnaroundTimeThreshold": 8
         |   }
         | }
         |}
         |""".stripMargin.parseJson

    val extras = Extras.parse(runSpec)
    extras.defaultTaskDxAttributes.get.treeTurnaroundTimeThreshold should be(Some(1))
    extras.perTaskDxAttributes.get.map {
      case (k, v) => k -> v.treeTurnaroundTimeThreshold
    } should be(
        Map(
            "Add" -> Some(2),
            "Multiply" -> Some(3)
        )
    )
    extras.defaultWorkflowDxAttributes.get.treeTurnaroundTimeThreshold should be(Some(5))
    extras.perWorkflowDxAttributes.get.map {
      case (k, v) => k -> v.treeTurnaroundTimeThreshold
    } should be(
        Map(
            "NestedWf" -> Some(8)
        )
    )
  }

  it should "include optional details in per task attributes" in {
    val runSpec: JsValue =
      """|{
         | "defaultTaskDxAttributes" : {
         |   "runSpec": {
         |     "timeoutPolicy": {
         |        "*": {
         |          "hours": 12
         |        }
         |     }
         |   }
         |  },
         | "perTaskDxAttributes" : {
         |   "Add": {
         |      "details": {
         |        "upstreamProjects": [
         |          {
         |            "name": "GATK4",
         |            "repoUrl": "https://github.com/broadinstitute/gatk",
         |            "version": "GATK-4.0.1.2",
         |            "license": "BSD-3-Clause",
         |            "licenseUrl": "https://github.com/broadinstitute/LICENSE.TXT",
         |            "author": "Broad Institute"
         |          }
         |        ]
         |      }
         |    }
         |  }
         |}
         |""".stripMargin.parseJson

    val extras = Extras.parse(runSpec)
    extras.defaultTaskDxAttributes should be(
        Some(
            DxAppJson(
                Some(DxRunSpec(None, None, None, Some(DxTimeout(None, Some(12), None)))),
                None
            )
        )
    )
    extras.perTaskDxAttributes should be(
        Some(
            Map(
                "Add" -> DxAppJson(
                    None,
                    Some(
                        DxDetails(
                            Some(
                                Vector(
                                    DxLicense("GATK4",
                                              "https://github.com/broadinstitute/gatk",
                                              "GATK-4.0.1.2",
                                              "BSD-3-Clause",
                                              "https://github.com/broadinstitute/LICENSE.TXT",
                                              "Broad Institute")
                                )
                            )
                        )
                    )
                )
            )
        )
    )
  }

  it should "parse the docker registry section" in {
    val data =
      """|{
         | "dockerRegistry" : {
         |   "registry" : "foo.bar.dnanexus.com",
         |   "username" : "perkins",
         |   "credentials" : "The Bandersnatch has gotten loose"
         | }
         |}
         |""".stripMargin.parseJson

    val extras = Extras.parse(data)
    extras.dockerRegistry should be(
        Some(
            DockerRegistry("foo.bar.dnanexus.com",
                           "The Bandersnatch has gotten loose",
                           Some("perkins"),
                           None)
        )
    )
  }

  it should "recognize errors in docker registry section" in {
    val data =
      """|{
         | "dockerRegistry" : {
         |   "registry_my" : "foo.bar.dnanexus.com",
         |   "username" : "perkins",
         |   "credentials" : "BandersnatchOnTheLoose"
         | }
         |}
         |""".stripMargin.parseJson
    assertThrows[DeserializationException] {
      Extras.parse(data)
    }
  }

  it should "recognize errors in docker registry section II" in {
    val data =
      """|{
         | "dockerRegistry" : {
         |   "registry" : "foo.bar.dnanexus.com",
         |   "credentials" : "BandersnatchOnTheLoose"
         | }
         |}
         |""".stripMargin.parseJson
    assertThrows[DeserializationException] {
      Extras.parse(data)
    }
  }

  it should "recognize errors in docker registry section III" in {
    val data =
      """|{
         | "dockerRegistry" : {
         |   "creds" : "XXX"
         | }
         |}
         |""".stripMargin.parseJson
    assertThrows[DeserializationException] {
      Extras.parse(data)
    }
  }

  it should "convert DxLicense to JsValue" in {
    val dxDetailsJson: JsValue =
      """|[
         |   {
         |      "author":"author1",
         |      "license":"license1",
         |      "licenseUrl":"licenseURL",
         |      "name":"name",
         |      "repoUrl":"repoURL",
         |      "version":"version1"
         |   },
         |   {
         |      "author":"author2",
         |      "license":"license2",
         |      "licenseUrl":"licenseURL",
         |      "name":"name2",
         |      "repoUrl":"repoURL",
         |      "version":"version2"
         |   }
         |]
         |""".stripMargin.parseJson

    val upstreamProjects = Vector(
        DxLicense("name", "repoURL", "version1", "license1", "licenseURL", "author1"),
        DxLicense("name2", "repoURL", "version2", "license2", "licenseURL", "author2")
    )

    val dxDetails = DxDetails(Some(upstreamProjects))

    val result = dxDetails.toJson
    result.asJsObject.fields("upstreamProjects") should be(dxDetailsJson)
  }

  it should "all DxAttr to return RunSpec Json" in {
    val expectedPolicy = """
                           |{
                           |  "*": {
                           |    "minutes": 30
                           |  }
                           |}
          """.stripMargin.parseJson

    val expected: Map[String, JsValue] = Map("timeoutPolicy" -> expectedPolicy)

    val dxAppJson: DxAppJson = DxAppJson(
        Some(DxRunSpec(None, None, None, Some(DxTimeout(None, None, Some(30))))),
        None
    )

    val runSpecJson: Map[String, JsValue] = dxAppJson.getApiRunSpecJson
    runSpecJson should be(expected)

  }

  it should "all DxAttr to return empty runSpec and details Json" in {
    val dxAppJson = DxAppJson(None, None)
    val runSpecJson = dxAppJson.getApiRunSpecJson
    runSpecJson should be(Map.empty)
    val detailJson = dxAppJson.getDetailsJson
    detailJson should be(Map.empty)
  }

  it should "all DxAttr to return Details Json" in {
    val expectedContent =
      """
        |[
        |  {
        |    "name": "GATK4",
        |    "repoUrl": "https://github.com/broadinstitute/gatk",
        |    "version": "GATK-4.0.1.2",
        |    "license": "BSD-3-Clause",
        |    "licenseUrl": "https://github.com/broadinstitute/LICENSE.TXT",
        |    "author": "Broad Institute"
        |  }
        |]
        |
          """.stripMargin.parseJson

    val expected: Map[String, JsValue] = Map("upstreamProjects" -> expectedContent)

    val dxAppJson: DxAppJson = DxAppJson(
        None,
        Some(
            DxDetails(
                Some(
                    Vector(
                        DxLicense("GATK4",
                                  "https://github.com/broadinstitute/gatk",
                                  "GATK-4.0.1.2",
                                  "BSD-3-Clause",
                                  "https://github.com/broadinstitute/LICENSE.TXT",
                                  "Broad Institute")
                    )
                )
            )
        )
    )

    val detailsJson = dxAppJson.getDetailsJson

    expected should be(detailsJson)
  }
}
