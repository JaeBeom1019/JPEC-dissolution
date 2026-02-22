# JPEC ForceFreeStates 변수/저장 경로 정리

실행 검증:
- 실제 실행: `JPEC.main(["/Users/jaebeomcho/Desktop/code/JPEC-dissolution/examples/DIIID-like_ideal_example"])`
- 출력/동작이 질문에서 주신 로그와 동일하게 재현됨
- 출력 파일 확인: `examples/DIIID-like_ideal_example/jpec.h5`

## 1) LCFS(마지막 경계)에서 Jacobian은 어디에 저장되나?

핵심 저장 위치는 **Equilibrium 객체의 `rzphi_jac`** 입니다.
- 코드: `equil.rzphi_jac` (`src/Equilibrium/EquilibriumTypes.jl`, `src/Equilibrium/DirectEquilibrium.jl`)
- HDF5 저장: `splines/rzphi/jac`

관련 파생 저장:
- ForceFreeStates metric 생성 시 `MetricData.fs[:,:,7] = J`, `fs[:,:,8] = dJ/dpsi` 로 복사됨 (`src/ForceFreeStates/Fourfit.jl`)
- 이후 Fourier 계수로 `jmat`, `jmat1` 사용 (`src/ForceFreeStates/Fourfit.jl`)

LCFS 값 접근 관점:
- equilibrium 그리드 마지막 면(보통 `psi ~ 1`)의 J는 `equil.rzphi_jac.nodal_derivs.partials[1, end, :]`
- ForceFreeStates 적분 경계(실제 `psilim`)는 `intr.psilim`이며, response 계산에서 boundary index는 `ForceFreeStates_results.step`를 사용

### 1-1) `nodal_derivs.partials`의 첫/둘째/셋째 엔트리 의미 (중요)

이건 코드에서 명시적으로 다음 의미로 사용됩니다 (`src/ForceFreeStates/Fourfit.jl` 주석 + 사용 패턴):
- `partials[1, :, :]` = 함수값 `f(psi, theta)`
- `partials[2, :, :]` = `∂f/∂psi`
- `partials[3, :, :]` = `∂f/∂theta`

즉 `equil.rzphi_jac.nodal_derivs.partials[1, ipsi, itheta]`에서:
- 첫 번째 index(`1`)는 "함수값"
- 두 번째 index(`ipsi`)는 `psi` grid index
- 세 번째 index(`itheta`)는 `theta` grid index

추가로, 실제 grid 자체는:
- `equil.rzphi_xs` = `psi` grid
- `equil.rzphi_ys` = normalized `theta` grid (`0..1`, 실제 각은 `2π*theta`)

### 1-2) `MetricData.fs[:,:,7]`의 축은 무엇인가?

`MetricData.fs` shape은 `(mpsi, mtheta, 8)`이고,
- 첫 번째 축 = `psi index`
- 두 번째 축 = `theta index`
- 세 번째 축 = 물리량 종류(1~8)

따라서 `MetricData.fs[:,:,7]`은 **`J(psi,theta)`** 입니다.  
`theta,zeta`가 아닙니다.

### 1-3) 그러면 `J`는 zeta independent 인가?

현재 JPEC equilibrium은 축대칭(2D) 표현이므로, 저장되는 Jacobian은 `J(psi,theta)`이며 **명시적 `zeta` 축이 없습니다**.

즉 코드 레벨에서는 `J = J(psi,theta)`로 취급하고, `zeta` 적분은 Fourier 직교성/주기성으로 분리됩니다.  
`papers/GPEC_notes.pdf` D절에서도 Fourier 전개 후 poloidal mode coupling에서 `theta` 적분 형태가 직접 등장합니다.

### 1-4) 이번 DIII-D run(`jac_type=hamada`)에서의 추가 사실

이 케이스는 `jac_type="hamada"`라 내부 지수들이 `(power_bp, power_b, power_r)=(0,0,0)`로 강제됩니다.  
실제 결과 파일에서 `splines/rzphi/jac`의 각 `psi` 행에 대해 `theta` 방향 표준편차를 확인하면:

- `max std over theta per psi row = 2.42e-13`

즉 수치적으로 이번 케이스는 사실상 `J(psi)`입니다 (theta 의존 거의 없음).  
하지만 일반적인 설정/좌표형(`jac_type`)에서는 `J(psi,theta)`가 될 수 있습니다.

## 2) 여기 pressure는 순수 pressure인가, normalization 된 값인가?

JPEC 내부 profile pressure는 **순수 P가 아니라 `mu0 * P`** 입니다.
- 정의: `P_spline: μ₀*P` (`src/Equilibrium/EquilibriumTypes.jl`)
- EFIT 읽기에서도 `pres_data .* mu0`로 변환 (`src/Equilibrium/ReadEquilibrium.jl`)
- HDF5 저장 키도 명시적으로 `splines/profiles/mu0p`

즉, 코드에서 `profiles.P_spline`/`profiles.P_deriv`는 기본적으로 `μ0P` 기준입니다.

## 3) psi/theta/zeta에 해당하는 Jacobian은 어떻게 저장되나?

이 코드에서 Jacobian은 **하나의 스칼라 장 `J(psi, theta)`** 로 저장됩니다.
- 저장: `equil.rzphi_jac` (`psi,theta` 2D)
- HDF5: `splines/rzphi/jac` (shape 예시: `(129,257)`)

중요:
- equilibrium은 축대칭 표현이라 Jacobian에 **명시적 zeta dependence를 따로 저장하지 않음**
- 따라서 `J_psi`, `J_theta`, `J_zeta`처럼 좌표별로 분리된 Jacobian 변수는 없음
- 대신 metric 조합으로 `g^ij * J` 계열을 저장/사용 (`MetricData.fs[:,:,1:6]`)
- `dJ/dpsi`는 `MetricData.fs[:,:,8]`에 저장됨

## 4) `dp/dpsi`는 어디에 저장? normalization 되었나?

저장/사용 위치:
- 기본: `equil.profiles.P_deriv` (spline derivative view)
- 사용 예: `pressure_gradient = profiles.P_deriv(psi)` in `src/ForceFreeStates/Bal.jl`

정규화 해석:
- `P` 자체가 `μ0P`
- 독립변수 `psi`는 normalized flux (`psi_norm`)
- 따라서 `profiles.P_deriv(psi)`는 **`d(μ0P)/d(psi_norm)`**

물리 flux `Ψ`에 대한 도함수가 필요하면:
- `d(μ0P)/dΨ = (d(μ0P)/dpsi_norm) / psio`
- `dP/dΨ = (d(μ0P)/dpsi_norm) / (μ0 * psio)`

참고:
- `dp/dpsi`를 별도 dataset으로 직접 쓰지는 않음(필요시 `mu0p` spline에서 미분해 복원)

## 5) boundary에서 각 theta, zeta에 대한 radial perturbation `xi_psi`는 어디에 저장?

두 가지 레벨이 있습니다.

1. ForceFreeStates 원본 eigenmode 저장:
- 변수: `odet.u_store[:,:,1,:]` (3번째 축 index 1이 `xi_psi`)
- HDF5: `integration/xi_psi` (shape 예시: `(34,34,136)`)
- boundary slice는 보통 `integration/xi_psi[:,:,end]` 또는 코드상 `ForceFreeStates_results.u_store[:,:,1,ForceFreeStates_results.step]`

2. Perturbed equilibrium 응답(가중합 후):
- HDF5: `perturbed_equilibrium/response/xi_psi_real`, `.../xi_psi_imag`
- shape 예시: `(npsi, mpert) = (136,34)`
- 이 값도 **(psi,m) 모드공간** 데이터

중요:
- 질문하신 "각 theta,zeta"의 실공간 `xi_psi(psi,theta,zeta)`를 직접 저장하는 경로는 현재 구현에 없음
- 현재는 Fourier mode coefficient 형태(`m,n`)로 저장하고, 필요시 후처리 inverse Fourier transform으로 `(theta,zeta)`를 복원해야 함

## 5-1) `Xi_psi`, `Xi_s`, `Xi_psi'`의 코드 내 실제 대응 (D절 식과 직접 매핑)

`papers/GPEC_notes.pdf` D절 식(49):
- `Xi_s = -A^{-1}(B Xi_psi' + C Xi_psi)`

코드(`src/ForceFreeStates/Sing.jl::sing_der!`) 대응:
- `u[:,:,1]` = `Xi_psi`
- `ud[:,:,1]` = `Xi_psi'`
- `ud[:,:,2]` = `Xi_s`

그리고 실제 계산은:
- `bmat <- A^{-1}B` (cholesky + `ldiv!`)
- `cmat <- A^{-1}C`
- `ud[:,:,2] = - bmat*ud[:,:,1] - cmat*u[:,:,1]`
- 즉 코드 구현은 정확히 `Xi_s = -A^{-1}(B Xi_psi' + C Xi_psi)`

주의:
- `u[:,:,2]`는 `Xi_s`가 아니라 보조 상태(코드 주석상 conjugate momenta 계열)입니다.
- HDF5에서:
  - `integration/xi_psi` = `u_store[:,:,1,:]`
  - `integration/dxi_psi` = `ud_store[:,:,1,:]`
  - `integration/xi_s` = `ud_store[:,:,2,:]`
  - `integration/u2` = `u_store[:,:,2,:]` (별도 상태)

## 5-2) `xi^s(psi,theta,zeta)`, `xi^psi(psi,theta,zeta)`를 edge에서 만들 수 있나?

결론: **가능합니다. 다만 현재 저장은 mode-space 계수이므로, 실공간 변환을 후처리로 해야 합니다.**

필요 데이터:
- edge `psi`에서 계수 벡터
  - `v_{m,n}`: `Xi_psi` 계수 (`integration/xi_psi[:,:,edge_index]` 또는 response 후 `xi_modes.psi[end,:]`)
  - `u_{m,n}`: `Xi_s` 계수 (`integration/xi_s[:,:,edge_index]`)
- 모드 번호 매핑: `info/mn_index`

실공간 복원:
- `xi(psi,theta,zeta) = Σ_m Σ_n coeff_{m,n}(psi) exp(i(m theta - n zeta))`
- 현재 예제는 `npert=1`이므로 사실상 `Σ_m`만 수행하고 `exp(-i n zeta)` 위상 곱 추가

즉, edge에서 `theta,zeta` 격자 데이터로 나타내는 것은 가능하지만, 코드가 기본으로 저장하는 포맷은 `(m,n)` Fourier 계수입니다.

## 8) edge에서 `|B|` 또는 `B^2`는 어디에 저장되나?

in-memory:
- `equil.eqfun_B(psi,theta)`가 `|B|`입니다 (`src/Equilibrium/EquilibriumTypes.jl`, `DirectEquilibrium.jl`)
- edge면은 `psi = equil.rzphi_xs[end]`에서 평가

HDF5 기본 저장:
- 현재 `JPEC.jl` 기본 출력은 `eqfun_B`를 파일에 쓰지 않습니다.
- 즉 `jpec.h5`에는 기본적으로 `splines/eqfun_B`가 없습니다.

따라서 edge `|B|`, `B^2`가 필요하면:
1. 실행 중 객체에서 `equil.eqfun_B((psi_edge, theta))`로 샘플링
2. 또는 코드에 dataset write를 추가해 저장

`B^2`는 단순히 `equil.eqfun_B(...)^2`로 얻을 수 있습니다.

## 9) `DivergenceFreeXi.tex` 기준: edge에서 `xi_theta`, `xi_zeta`에 필요한 입력값 정리

`papers/DivergenceFreeXi.tex` 식 기준으로 edge에서 필요한 값은 다음입니다.

필수 입력(모드공간):
- `Xi_psi` at edge: `odet.u_store[:, solution_index, 1, edge_step]`
- `Xi_psi'` at edge: `odet.ud_store[:, solution_index, 1, edge_step]`
- `A,B,C` at edge psi: `ffit.amats/bmats/cmats`를 `psi_edge`에서 평가
- `m` 범위: `intr.mlow:intr.mhigh`
- `n` 모드: `intr.nlow:intr.nhigh` (현재 예제는 single-`n`)

필수 입력(실공간):
- `theta` grid: `equil.rzphi_ys` (normalized `[0,1]`)
- edge Jacobian: `J(psi_edge,theta)` from `equil.rzphi_jac((psi_edge,theta))`
- 이전 psi 면 Jacobian/`xi_psi`:
  - `psi_prev = odet.psi_store[edge_step-1]`
  - `J(psi_prev,theta)`, `xi_psi(psi_prev,theta)`
  - 이유: `∂ψ(J xi^psi)` 계산에 필요
- edge `q`: `q_edge = odet.q_store[edge_step]`
- `chi'`: 코드에서 `chi1 = 2π * equil.psio`

정의/식 매핑:
- `Xi_s = -A^{-1}(B Xi_psi' + C Xi_psi)` (`tex` Eq.49)
- `S = -∂ψ(J xi^psi) - (2π i n / chi') J xi^s`
  - normalized angle에서 phase를 `exp(2π i (mθ-nζ))`로 두면 `2π`가 붙음
- `f = J xi^theta`에 대해 `f' - (2π i n q) f = S`
- `xi^zeta = q xi^theta - xi^s / chi'`

## 10) 기존 문서에 없었던 항목 보강

이 문서에 추가된 핵심 누락 항목:
- `A,B,C`를 edge에서 재평가해 `Xi_s`를 직접 재구성하는 경로
- `xi_theta` 계산에서 필요한 `∂ψ(J xi^psi)` 때문에 edge 한 점만으로는 부족하고, 직전 psi 저장점이 반드시 필요하다는 점
- `theta` grid가 `equil.rzphi_ys`라는 점(0~1 normalized)
- `chi' = 2π*psio`의 코드 대응
- `exp(2π i(mθ-nζ))` 기준으로 real-space 재구성/적분해야 한다는 구현상 기준

## 11) 구현 방향 및 현재 반영 상태 (`Surface.jl`)

요청대로 기존 핵심 파일(`EulerLagrange.jl`, `Fourfit.jl` 등)은 수정하지 않고:
- 신규 파일 추가: `src/ForceFreeStates/Surface.jl`
- include 추가: `src/ForceFreeStates/ForceFreeStates.jl`에 `include("Surface.jl")`

`Surface.jl` 제공 함수:
- `compute_edge_xi_components(equil, ffit, odet, intr; solution_index=1, n_mode=intr.nlow, zeta=0.0)`
  - edge에서 `Xi_psi`, `Xi_psi'`, `A,B,C`로 `Xi_s` 계산
  - `theta`별 `xi_psi`, `xi_s`, `xi_theta`, `xi_zeta` 복원
  - `Xi_s(ABC)`와 저장된 `Xi_s(ud_store[:,:,2])` 차이도 반환
- `print_edge_xi_components(...)`
  - 각 theta 위치에서 `xi_psi`, `xi_theta`, `xi_zeta`를 출력

구현 세부 방향(반영됨):
1. edge/직전 psi를 사용한 backward difference로 `∂ψ(J xi^psi)` 계산
2. `theta ∈ [0,1]` grid에서 복소 Trapezoidal 적분으로 ODE 직접 풂
3. periodic 조건으로 integration constant(`f(0)`) 결정
4. 마지막에 `xi_zeta = q xi_theta - xi_s/chi'` 계산

추가 주의점:
- `e^{inq}=1` (rational condition 근처)에서는 분모 `1-e^{inq}`가 작아질 수 있어 경고 필요
- 현재 예제(`npert=1`)는 단일 `n`; multi-`n`에서는 원하는 `n_mode` block을 선택해서 동일 절차 적용

동작 검증:
- 실제 run에서 `compute_edge_xi_components(...)` 실행 확인
- `max |Xi_s(ABC) - Xi_s(stored)| ≈ 1.0e-16` 수준으로 일치 확인

## 6) ForceFreeStates 디렉토리 파일 구조(질문 항목과 직접 관련된 부분)

- `src/ForceFreeStates/ForceFreeStatesStructs.jl`
  - 핵심 상태 구조체: `ForceFreeStatesInternal`, `FourFitVars`, `OdeState`
  - `OdeState.u_store/ud_store`에 적분 결과 저장
- `src/ForceFreeStates/Fourfit.jl`
  - metric/Jacobian(`J`, `dJ/dpsi`) 계산 및 Fourier-matrix 구성
- `src/ForceFreeStates/EulerLagrange.jl`
  - ODE 적분, 저장 step 관리, singular crossing
- `src/ForceFreeStates/Bal.jl`
  - `pressure_gradient = profiles.P_deriv(psi)` 사용
- `src/ForceFreeStates/FixedBoundaryStability.jl`
  - 저장된 적분해로 fixed-boundary criterion 평가
- `src/ForceFreeStates/Free.jl`
  - vacuum/free-boundary 에너지 처리

## 7) 확인한 실제 HDF5 shape (DIIID example)

- `integration/xi_psi`: `(34, 34, 136)`
- `integration/psi`: `(136,)`
- `splines/rzphi/jac`: `(129, 257)`
- `splines/profiles/mu0p`: `(129,)`
- `perturbed_equilibrium/response/xi_psi_real`: `(136, 34)`

## 12) `DivergenceFreeXi.tex` 기준: `∇·xi`를 계산할 때 필요한 것

대상 식:
- `∇·xi = (1/J) [ ∂ψ(J xi^ψ) + ∂θ(J xi^θ) + ∂ζ(J xi^ζ) ]`
- 단일 `n`의 normalized 각도 정의(`exp(2π i(mθ-nζ))`)를 쓰면 `∂ζ -> -2π i n`

필요 입력값:
- `xi^ψ(psi,theta,zeta)` 또는 mode-space `Xi_psi` + inverse Fourier
- `xi^θ(psi,theta,zeta)` (현재 `Surface.jl`에서 복원)
- `xi^ζ(psi,theta,zeta)` (현재 `Surface.jl`에서 복원)
- `J(psi,theta)` (`equil.rzphi_jac`)
- `∂ψ(J xi^ψ)` 계산용 `psi` 인접점 데이터 (`edge`만 볼 때는 `psi_edge`, `psi_prev`)
- `q(psi)` (`equil.profiles.q_spline` 또는 `odet.q_store`)
- `n` 모드번호 (`intr.nlow:intr.nhigh`)
- `chi' = 2π*psio` (`chi1 = 2π*equil.psio`)
- `theta` grid (`equil.rzphi_ys`)

코드에서 바로 확보 가능한 것:
- `Xi_psi`: `odet.u_store[:,:,1,:]`
- `Xi_psi'`: `odet.ud_store[:,:,1,:]`
- `Xi_s`: `odet.ud_store[:,:,2,:]`
- `J`: `equil.rzphi_jac((psi,theta))`
- `q`: `odet.q_store[...]` 또는 `equil.profiles.q_spline(psi)`
- `theta` grid: `equil.rzphi_ys`

실무 포인트:
- `∇·xi=0`는 현재 `xi_theta` 복원 ODE 구성에 이미 반영되어 있어 해석적으로는 항등적으로 0이 되도록 설계됨 (`DivergenceFreeXi.tex` Step 4).
- 실제 수치 검증에서는 차분 오차 때문에 작은 residual이 생길 수 있으므로, 절대값과 상대값 둘 다 보는 게 안전함.

## 13) `DivergenceFreeXi.tex` 기준: `xi_perp`를 계산할 때 필요한 것

대상 식:
- `xi_perp = xi - ((xi·B)/B^2) B`
- 성분식으로는 `xi_parallel = (xi·B)/B^2`를 먼저 구한 뒤 `xi^θ, xi^ζ`에서 `B^θ, B^ζ` 방향 성분을 제거

필요 입력값:
- `xi^ψ, xi^θ, xi^ζ` (real-space)
- `B`의 contravariant 성분:
  - `B^ψ = 0`
  - `B^θ = chi'/J`
  - `B^ζ = q*chi'/J`
- `B^2`
- `xi·B`
- 공변 metric tensor `g_ij` (또는 동등한 기저벡터 정보)

코드에서 바로 있는 것:
- `J`, `q`, `chi'`는 이미 확보 가능
- `|B|`는 `equil.eqfun_B(psi,theta)`로 확보 가능 (`B^2 = eqfun_B^2`)
- `xi^ψ, xi^θ, xi^ζ`는 `Surface.jl` 복원 결과로 확보 가능

현재 바로 안 나오는 것(중요):
- `g_ij` 전체는 HDF5에 직접 저장되어 있지 않음
- `eqfun_metric1`, `eqfun_metric2`만으로는 `g_ij` 전체를 단독 복원하기 부족

`g_ij`를 얻는 현실적 방법:
- 방법 A: `rzphi` 미분값(`partials`)으로 basis를 만들어 현장에서 `g^ij` 계산 후 3x3 inverse로 `g_ij` 복원
  - `Fourfit.jl`, `Bal.jl`의 `v`/`w` 계산식 재사용 가능
- 방법 B: `MetricData.fs[:,:,1:6] = g^{ij}*J`를 재평가/재구성해 `g^{ij}`를 만들고 inverse

`xi_perp` 계산 시 필요한 최소 개발 항목(미구현):
- `Surface.jl`에서 `(psi,theta)`별 `g^{ij}` 평가 함수
- 3x3 역행렬로 `g_ij` 구성
- `xi·B`, `B^2`, `xi_parallel`, `xi_perp` 성분/크기 계산 함수
- 필요 시 edge 전용 출력/플롯(`|xi_perp|`, `xi_parallel`) 함수

## 14) `DivergenceFreeXi.tex` Step 9 기준 `deltaWs_full`에 필요한 항목 (정리 + 코드 매핑)

목표 식:
- `F = div(xi_perp) + 2 xi_perp·kappa`
- `deltaWs_full = -1/2 ∫_0^1 dtheta J xi_psi^* B^2 F`  (single-`n`에서 `∫dζ=1`)

필요 물리량과 코드 소스:
- `xi^psi, xi^theta, xi^zeta`:
  - `Surface.jl::compute_edge_xi_components(...)` 출력 사용
- `J(psi,theta)`:
  - `data.J_edge` (`equil.rzphi_jac` 평가)
- `q(psi_edge)`:
  - `data.q_edge` (`odet.q_store` 기반)
- `chi'`:
  - `chi1 = 2π * equil.psio`
- `p'`:
  - `pprime_mu0 = equil.profiles.P_deriv(psi_edge)` = `d(mu0*P)/d(psi_norm)`
- `B^2`, `∂psi B^2`, `∂theta B^2`:
  - `equil.eqfun_B((psi,theta))`, `deriv=(1,0)`, `deriv=(0,1)`로 계산
  - `B^2 = B^2`, `dB2_dpsi = 2B*dB_dpsi`, `dB2_dtheta = 2B*dB_dtheta`
- metric tensor:
  - `rzphi_*` partials에서 `v`(contravariant basis) 구성
  - `g^{ij}` 계산 후 `inv(g^{ij}) -> g_ij`

구현 위치:
- `src/ForceFreeStates/Surface.jl::compute_edge_surface_energy_terms(...)`

## 15) `deltaWs_full` 구현 알고리즘 (`Surface.jl`)과 수식 일치성 체크

실제 구현 순서:
1. edge에서 `xi^psi, xi^theta, xi^zeta` 복원 (`compute_edge_xi_components`)
2. 각 `theta`에서 geometry 재구성:
   - `rzphi_rsquared/offset/nu/jac`와 `psi/theta` 미분값으로 `g^{ij}` 생성
   - `g_ij = inv(g^{ij})`
3. `B`와 도함수 평가:
   - `B^2`는 Bal/Fourfit 일관식 사용:
     - `B^2 = chi'^2 * (v21^2 + v22^2 + (v23 + q*v33)^2)`
   - 도함수는 `eqfun_B` spline 도함수로 계산:
     - `dB2_dpsi = 2B*dB_dpsi`, `dB2_dtheta = 2B*dB_dtheta`
4. `xi_parallel/B` 계산:
   - `lambda = (xi·B)/B^2`
5. `xi_perp` 계산:
   - `xi_perp = xi - lambda * B`
6. `div(xi_perp)` 계산 (normalized theta):
   - `div_xi_perp = -(chi'/J) * (∂theta - 2π i n q) (xi_parallel/B)`
   - `∂theta`는 periodic central difference
7. Step 7b direct 식으로 `xi_perp·kappa` 계산:
   - `xi_perp·kappa = [ (pprime_mu0 + 0.5*dB2_dpsi)*(xi_perp·∇psi) + 0.5*dB2_dtheta*(xi_perp·∇theta) ] / B^2`
8. `F = div_xi_perp + 2*xi_perp·kappa`
9. 적분:
   - `es_discontinuous = -(1/(2*mu0)) * pprime_mu0 * ∫ J|xi_psi|^2 dtheta`
   - `es_full = -0.5 * ∫ J * conj(xi_psi) * B^2 * F dtheta`

수식-코드 consistency 체크 포인트:
- `pprime_mu0 = d(mu0*P)/dpsi`를 그대로 사용하고, `es_full`은 사용자 정의에 맞게 `B^2` 항만 포함
- 단일 `n` 모드 기준으로 `∫dζ=1` 처리
- `theta`는 normalized grid (`equil.rzphi_ys`), 따라서 toroidal/poloidal phase operator에 `2π`가 포함됨

최종 출력 함수:
- `Surface.jl::print_edge_pressure_drive_comparison(...)`
  - 현재는 `ep, et, es_discontinuous, es_full`을 모두 출력
  - 실수부 비교(`Re(es_*) / Re(ep or et)`)도 함께 출력
- alias:
  - `print_edge_surface_energy_comparison(...)`
