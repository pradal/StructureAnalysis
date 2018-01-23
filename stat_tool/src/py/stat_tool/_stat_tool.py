__all__ = []


# Import Boost.Python module
from . import __stat_tool

# Resolve scopes
__stat_tool.std.IosBase.event = __stat_tool.std._ios_base.event
__stat_tool.std.IosBase.Init = __stat_tool.std._ios_base.Init
__stat_tool.std.Locale.Facet = __stat_tool.std._locale.Facet
__stat_tool.std._BasicOstream_e1391944268253558f04b6f996bb5a8b.Sentry = __stat_tool.std.__basic_ostream_e1391944268253558f04b6f996bb5a8b.Sentry
__stat_tool.std.Locale.Id = __stat_tool.std._locale.Id
__stat_tool.std.IosBase.Failure = __stat_tool.std._ios_base.Failure

# Group template specializations
__stat_tool.std._Allocator = (__stat_tool.std._Allocator_0e202a8855575dee86e8a8c24e73a6c0, __stat_tool.std._Allocator_1663f1821363576d83bbdfbc545ed162, __stat_tool.std._Allocator_1f29b58eee9156678af1ee8254be47f4, __stat_tool.std._Allocator_8180e493d9f4506e83585d89b1503ab0, __stat_tool.std._Allocator_a8a9526a25e8542295dd2df108ca409c, __stat_tool.std._Allocator_bc5c55c7605154049e618eb1601ff1eb, __stat_tool.std._Allocator_e01b7ade8cab5e31a34a5e2f80a619a6, __stat_tool.std._Allocator_f2215828ca3a504b8e3d080988e9947f)
__stat_tool.std._BasicStreambuf = (__stat_tool.std._BasicStreambuf_112dc12b863f53fea4df7b3ba388fd84)
__stat_tool.std._CharTraits = (__stat_tool.std._CharTraits_277a0516fe4451448165550d8b9d6b2b)
__stat_tool.std._BasicIos = (__stat_tool.std._BasicIos_f8b8546034205658b6e3e16175284f26)
__stat_tool.stat_tool._TemplateMultiPlotSet = (__stat_tool.stat_tool._TemplateMultiPlotSet_4f445857b49f50029570f15c5ecf9324)
__stat_tool.std._Ctype = (__stat_tool.std._Ctype_488de6b23c2d582c8382ac19e518b6a8)
__stat_tool.std._InitializerList = (__stat_tool.std._InitializerList_0155e43161935b21b27e7d2d4f43340d, __stat_tool.std._InitializerList_46156c79cdee52329533e5a44157e2e1, __stat_tool.std._InitializerList_4b65c81e2b1257af808046c0032f49ce, __stat_tool.std._InitializerList_55e5f2d64ff95ddc8d2f941abc315210, __stat_tool.std._InitializerList_821f52165ab1565093f9ea5b7cdd5401, __stat_tool.std._InitializerList_b6671fe8acff574aaf15e1ac26abb1bf, __stat_tool.std._InitializerList_bd3aca0a6cb156939d673e0f13c9dca5, __stat_tool.std._InitializerList_cc80f44cefab5e0f9c9fff48bf8aa217)
__stat_tool.std._MoveIterator = (__stat_tool.std._MoveIterator_f38d1184226e5eab899f2ead9c835dc3)
__stat_tool.std._BasicString = (__stat_tool.std._BasicString_448c20257e485acda59dc59305fceb58)
__stat_tool.std._Vector = (__stat_tool.std._Vector_107131f9768c56e794a9b0de728d1738, __stat_tool.std._Vector_6b9ae5eac40858c9a0f5e6e21c15d1d3, __stat_tool.std._Vector_81c18c08844552cc93122f9c501ecc0d, __stat_tool.std._Vector_c5260d1a7cbb5ef3bc32b582df09a801, __stat_tool.std._Vector_d3fbcd4a393754ca8d7be71823564225, __stat_tool.std._Vector_db7f50b235e15b1aa024911715fa604a, __stat_tool.std._Vector_f0840635c9e2594cb061b2ce1bc3514a)
__stat_tool.std._IntegralConstant = (__stat_tool.std._IntegralConstant_93823bb50db15ce2a0108011ea943902)
__stat_tool.std._Pair = (__stat_tool.std._Pair_15bfc45e4cb45511b08b9e71d0c7f189, __stat_tool.std._Pair_ad547ba62d555da5925f982f13f6787b)
__stat_tool.std._BasicOstream = (__stat_tool.std._BasicOstream_e1391944268253558f04b6f996bb5a8b)
__stat_tool.stat_tool._ChainReestimation = (__stat_tool.stat_tool._ChainReestimation_4a0641faa57256529b8c34bd13a4c984)
__stat_tool.stat_tool._Reestimation = (__stat_tool.stat_tool._Reestimation_4a7eb3f23eb959139a416b2a8b293302, __stat_tool.stat_tool._Reestimation_710d7ee5573c5d7f8f9127a08b4f3dfd)

# Define aliases
__stat_tool.stat_tool.Range = __stat_tool.std._Pair_ad547ba62d555da5925f982f13f6787b
__stat_tool.stat_tool.Label = __stat_tool.std._Pair_15bfc45e4cb45511b08b9e71d0c7f189
__stat_tool.stat_tool.MultiPlotSet = __stat_tool.stat_tool._TemplateMultiPlotSet_4f445857b49f50029570f15c5ecf9324
__stat_tool.stat_tool.PlotPoint = __stat_tool.std._Pair_ad547ba62d555da5925f982f13f6787b
__stat_tool.std._BasicStreambuf_112dc12b863f53fea4df7b3ba388fd84.TraitsType = __stat_tool.std._CharTraits_277a0516fe4451448165550d8b9d6b2b
__stat_tool.std._BasicStreambuf_112dc12b863f53fea4df7b3ba388fd84.StreambufType = __stat_tool.std._BasicStreambuf_112dc12b863f53fea4df7b3ba388fd84
__stat_tool.std._Pair_15bfc45e4cb45511b08b9e71d0c7f189.FirstType = __stat_tool.std._Pair_ad547ba62d555da5925f982f13f6787b
__stat_tool.std._Pair_15bfc45e4cb45511b08b9e71d0c7f189.SecondType = __stat_tool.std._BasicString_448c20257e485acda59dc59305fceb58
__stat_tool.std._BasicString_448c20257e485acda59dc59305fceb58.TraitsType = __stat_tool.std._CharTraits_277a0516fe4451448165550d8b9d6b2b
__stat_tool.std._InitializerList_46156c79cdee52329533e5a44157e2e1.ValueType = __stat_tool.stat_tool.discrete_parametric
__stat_tool.std._InitializerList_55e5f2d64ff95ddc8d2f941abc315210.ValueType = __stat_tool.stat_tool.FrequencyDistribution
__stat_tool.std.IosBase.Openmode = __stat_tool.std.ios_openmode
__stat_tool.std._IntegralConstant_93823bb50db15ce2a0108011ea943902.Type = __stat_tool.std._IntegralConstant_93823bb50db15ce2a0108011ea943902
__stat_tool.std._InitializerList_b6671fe8acff574aaf15e1ac26abb1bf.ValueType = __stat_tool.stat_tool.Vectors
__stat_tool.std._InitializerList_bd3aca0a6cb156939d673e0f13c9dca5.ValueType = __stat_tool.stat_tool.DiscreteParametric
__stat_tool.std._InitializerList_cc80f44cefab5e0f9c9fff48bf8aa217.ValueType = __stat_tool.stat_tool.process_distribution
__stat_tool.std._MoveIterator_f38d1184226e5eab899f2ead9c835dc3.IteratorCategory = __stat_tool.std.RandomAccessIteratorTag
__stat_tool.std._MoveIterator_f38d1184226e5eab899f2ead9c835dc3.ValueType = __stat_tool.stat_tool.discrete_parametric
__stat_tool.std._Allocator_0e202a8855575dee86e8a8c24e73a6c0.PropagateOnContainerMoveAssignment = __stat_tool.std._IntegralConstant_93823bb50db15ce2a0108011ea943902
__stat_tool.std._Vector_107131f9768c56e794a9b0de728d1738.AllocatorType = __stat_tool.std._Allocator_1663f1821363576d83bbdfbc545ed162
__stat_tool.std._Allocator_1663f1821363576d83bbdfbc545ed162.PropagateOnContainerMoveAssignment = __stat_tool.std._IntegralConstant_93823bb50db15ce2a0108011ea943902
__stat_tool.std._Allocator_1f29b58eee9156678af1ee8254be47f4.ValueType = __stat_tool.stat_tool.DiscreteParametric
__stat_tool.std._Allocator_1f29b58eee9156678af1ee8254be47f4.PropagateOnContainerMoveAssignment = __stat_tool.std._IntegralConstant_93823bb50db15ce2a0108011ea943902
__stat_tool.std._Vector_6b9ae5eac40858c9a0f5e6e21c15d1d3.AllocatorType = __stat_tool.std._Allocator_e01b7ade8cab5e31a34a5e2f80a619a6
__stat_tool.std._Allocator_8180e493d9f4506e83585d89b1503ab0.ValueType = __stat_tool.stat_tool.process_distribution
__stat_tool.std._Allocator_8180e493d9f4506e83585d89b1503ab0.PropagateOnContainerMoveAssignment = __stat_tool.std._IntegralConstant_93823bb50db15ce2a0108011ea943902
__stat_tool.std._Vector_81c18c08844552cc93122f9c501ecc0d.ValueType = __stat_tool.stat_tool.process_distribution
__stat_tool.std._Vector_81c18c08844552cc93122f9c501ecc0d.AllocatorType = __stat_tool.std._Allocator_8180e493d9f4506e83585d89b1503ab0
__stat_tool.std._Allocator_a8a9526a25e8542295dd2df108ca409c.ValueType = __stat_tool.stat_tool.discrete_parametric
__stat_tool.std._Allocator_a8a9526a25e8542295dd2df108ca409c.PropagateOnContainerMoveAssignment = __stat_tool.std._IntegralConstant_93823bb50db15ce2a0108011ea943902
__stat_tool.std._Allocator_bc5c55c7605154049e618eb1601ff1eb.ValueType = __stat_tool.stat_tool.Vectors
__stat_tool.std._Allocator_bc5c55c7605154049e618eb1601ff1eb.PropagateOnContainerMoveAssignment = __stat_tool.std._IntegralConstant_93823bb50db15ce2a0108011ea943902
__stat_tool.std._Vector_c5260d1a7cbb5ef3bc32b582df09a801.ValueType = __stat_tool.stat_tool.DiscreteParametric
__stat_tool.std._Vector_c5260d1a7cbb5ef3bc32b582df09a801.AllocatorType = __stat_tool.std._Allocator_1f29b58eee9156678af1ee8254be47f4
__stat_tool.std._Vector_d3fbcd4a393754ca8d7be71823564225.ValueType = __stat_tool.stat_tool.discrete_parametric
__stat_tool.std._Vector_d3fbcd4a393754ca8d7be71823564225.AllocatorType = __stat_tool.std._Allocator_a8a9526a25e8542295dd2df108ca409c
__stat_tool.std._Vector_db7f50b235e15b1aa024911715fa604a.ValueType = __stat_tool.stat_tool.Vectors
__stat_tool.std._Vector_db7f50b235e15b1aa024911715fa604a.AllocatorType = __stat_tool.std._Allocator_bc5c55c7605154049e618eb1601ff1eb
__stat_tool.std._Allocator_e01b7ade8cab5e31a34a5e2f80a619a6.PropagateOnContainerMoveAssignment = __stat_tool.std._IntegralConstant_93823bb50db15ce2a0108011ea943902
__stat_tool.std._Vector_f0840635c9e2594cb061b2ce1bc3514a.ValueType = __stat_tool.stat_tool.FrequencyDistribution
__stat_tool.std._Vector_f0840635c9e2594cb061b2ce1bc3514a.AllocatorType = __stat_tool.std._Allocator_f2215828ca3a504b8e3d080988e9947f
__stat_tool.std._Allocator_f2215828ca3a504b8e3d080988e9947f.ValueType = __stat_tool.stat_tool.FrequencyDistribution
__stat_tool.std._Allocator_f2215828ca3a504b8e3d080988e9947f.PropagateOnContainerMoveAssignment = __stat_tool.std._IntegralConstant_93823bb50db15ce2a0108011ea943902
__stat_tool.std._BasicIos_f8b8546034205658b6e3e16175284f26.TraitsType = __stat_tool.std._CharTraits_277a0516fe4451448165550d8b9d6b2b
__stat_tool.std._BasicIos_f8b8546034205658b6e3e16175284f26.CtypeType = __stat_tool.std._Ctype_488de6b23c2d582c8382ac19e518b6a8
__stat_tool.std._BasicOstream_e1391944268253558f04b6f996bb5a8b.TraitsType = __stat_tool.std._CharTraits_277a0516fe4451448165550d8b9d6b2b
__stat_tool.std._BasicOstream_e1391944268253558f04b6f996bb5a8b.StreambufType = __stat_tool.std._BasicStreambuf_112dc12b863f53fea4df7b3ba388fd84
__stat_tool.std._BasicOstream_e1391944268253558f04b6f996bb5a8b.IosType = __stat_tool.std._BasicIos_f8b8546034205658b6e3e16175284f26
__stat_tool.std._BasicOstream_e1391944268253558f04b6f996bb5a8b.OstreamType = __stat_tool.std._BasicOstream_e1391944268253558f04b6f996bb5a8b
__stat_tool.std._BasicOstream_e1391944268253558f04b6f996bb5a8b.CtypeType = __stat_tool.std._Ctype_488de6b23c2d582c8382ac19e518b6a8
