; ModuleID = 'probe4.4f73e7baea3a4be-cgu.0'
source_filename = "probe4.4f73e7baea3a4be-cgu.0"
target datalayout = "e-m:o-i64:64-i128:128-n32:64-S128"
target triple = "arm64-apple-macosx11.0.0"

@alloc_a4458fba36b2c9b9cbdb4250acf3d8dc = private unnamed_addr constant <{ [75 x i8] }> <{ [75 x i8] c"/rustc/190f4c96116a3b59b7de4881cfec544be0246d84/library/core/src/num/mod.rs" }>, align 1
@alloc_4ff0eb3999050a3ac98dc327b45b2d5d = private unnamed_addr constant <{ ptr, [16 x i8] }> <{ ptr @alloc_a4458fba36b2c9b9cbdb4250acf3d8dc, [16 x i8] c"K\00\00\00\00\00\00\00y\04\00\00\05\00\00\00" }>, align 8
@str.0 = internal unnamed_addr constant [25 x i8] c"attempt to divide by zero"

; probe4::probe
; Function Attrs: uwtable
define void @_ZN6probe45probe17had7f6015aa408dc0E() unnamed_addr #0 {
start:
  %0 = call i1 @llvm.expect.i1(i1 false, i1 false)
  br i1 %0, label %panic.i, label %"_ZN4core3num21_$LT$impl$u20$u32$GT$10div_euclid17h46a2ab4c060d8f3aE.exit"

panic.i:                                          ; preds = %start
; call core::panicking::panic
  call void @_ZN4core9panicking5panic17h10b6756e897165ceE(ptr align 1 @str.0, i64 25, ptr align 8 @alloc_4ff0eb3999050a3ac98dc327b45b2d5d) #3
  unreachable

"_ZN4core3num21_$LT$impl$u20$u32$GT$10div_euclid17h46a2ab4c060d8f3aE.exit": ; preds = %start
  ret void
}

; Function Attrs: nocallback nofree nosync nounwind willreturn memory(none)
declare i1 @llvm.expect.i1(i1, i1) #1

; core::panicking::panic
; Function Attrs: cold noinline noreturn uwtable
declare void @_ZN4core9panicking5panic17h10b6756e897165ceE(ptr align 1, i64, ptr align 8) unnamed_addr #2

attributes #0 = { uwtable "frame-pointer"="non-leaf" "probe-stack"="inline-asm" "target-cpu"="apple-m1" }
attributes #1 = { nocallback nofree nosync nounwind willreturn memory(none) }
attributes #2 = { cold noinline noreturn uwtable "frame-pointer"="non-leaf" "probe-stack"="inline-asm" "target-cpu"="apple-m1" }
attributes #3 = { noreturn }

!llvm.module.flags = !{!0}
!llvm.ident = !{!1}

!0 = !{i32 8, !"PIC Level", i32 2}
!1 = !{!"rustc version 1.77.0-nightly (190f4c961 2024-01-09)"}
