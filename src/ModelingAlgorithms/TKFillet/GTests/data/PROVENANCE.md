# Fillet GTest fixtures

These BREP fixtures are checked in so the tests do not depend on a local
labor directory. The paths and SHA-256 values below identify their sources.

| Fixture | Source | SHA-256 |
| --- | --- | --- |
| `PeriodicSeamSynthetic.brep` | Generated from OCCT primitives by [`make_PeriodicSeamSynthetic.tcl`](make_PeriodicSeamSynthetic.tcl); full R20 cylinder with a quarter-cylinder boss | `AF03872EB5537ABCE6EA1CC94D867CF15649E139FBAF70C10E4145ADDA80F1A5` |
| `bug28354.brep` | [`a-betenev/opencascade-dataset`, revision `d78a7a9e`](https://github.com/a-betenev/opencascade-dataset/blob/d78a7a9e6de57699261fa03641ec1b5bafac39a1/brep/bug28354.brep) | `847AB883D542419D01C68E3A45B482B5E3A75D6B2ECBBDD05E957AC15C781253` |
| `OCC602.brep` | [`a-betenev/opencascade-dataset`, revision `d78a7a9e`](https://github.com/a-betenev/opencascade-dataset/blob/d78a7a9e6de57699261fa03641ec1b5bafac39a1/brep/OCC602.brep) | `693F7A5D8A0A78871402788CFFD04B9D2D7D8F0B9A6E8A3CDC6D603388123ED7` |
| `bug26655_OCC50.brep` | [`a-betenev/opencascade-dataset`, revision `d78a7a9e`](https://github.com/a-betenev/opencascade-dataset/blob/d78a7a9e6de57699261fa03641ec1b5bafac39a1/brep/bug26655_OCC50.brep) | `C394355D83021B35FA874C397A028B55119BD7792BF22E93760B697ECFE1B8F0` |
| `Issue1761_YRib.brep` | Public [`Open-Cascade-SAS/OCCT` issue #1761](https://github.com/Open-Cascade-SAS/OCCT/issues/1761), attachment `conge-y-nervure-base.brep.txt` | `985CB51C2395CD7C67E670796A8FCE3411B092E5CF75E2492B1B7F983FB691C2` |

The dataset and issue attachments are third-party materials. Their source and
rights should be identified in any later upstream contribution under the
repository's contributor agreement. The synthetic seam fixture is regenerated
byte-for-byte by its Tcl script. Its Boolean topology differs from the separate
customer acceptance model, so the four public seam tests are independent
regressions rather than a replacement for that local acceptance run.
