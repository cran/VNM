# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

g1 <- function(T, order) {
    .Call(`_VNM_g1`, T, order)
}

g2 <- function(T, dt, order) {
    .Call(`_VNM_g2`, T, dt, order)
}

f234 <- function(T, x, order) {
    .Call(`_VNM_f234`, T, x, order)
}

SDM <- function(M) {
    .Call(`_VNM_SDM`, M)
}

Trans <- function(M) {
    .Call(`_VNM_Trans`, M)
}

Minus <- function(M1, M2) {
    .Call(`_VNM_Minus`, M1, M2)
}

Plus <- function(M1, M2) {
    .Call(`_VNM_Plus`, M1, M2)
}

Multiple <- function(M1, M2) {
    .Call(`_VNM_Multiple`, M1, M2)
}

infor234 <- function(T, x, order) {
    .Call(`_VNM_infor234`, T, x, order)
}

d1 <- function(T, x, xl, inv, order) {
    .Call(`_VNM_d1`, T, x, xl, inv, order)
}

d2 <- function(T, x, xl, inv, order) {
    .Call(`_VNM_d2`, T, x, xl, inv, order)
}

d3 <- function(T, x, xl, inv, dt, order) {
    .Call(`_VNM_d3`, T, x, xl, inv, dt, order)
}

dd1 <- function(T, x1, x2, xl, inv, order) {
    .Call(`_VNM_dd1`, T, x1, x2, xl, inv, order)
}

dd2 <- function(T, x1, x2, xl, inv, order) {
    .Call(`_VNM_dd2`, T, x1, x2, xl, inv, order)
}

dd3 <- function(T, x1, x2, xl, inv, dt, order) {
    .Call(`_VNM_dd3`, T, x1, x2, xl, inv, dt, order)
}

ds1 <- function(T, x, inv, order) {
    .Call(`_VNM_ds1`, T, x, inv, order)
}

ds2 <- function(T, x, inv, order) {
    .Call(`_VNM_ds2`, T, x, inv, order)
}

ds3 <- function(T, x, inv, dt, order) {
    .Call(`_VNM_ds3`, T, x, inv, dt, order)
}

sMultiple <- function(s, M) {
    .Call(`_VNM_sMultiple`, s, M)
}

upinfor <- function(W, T, X, order) {
    .Call(`_VNM_upinfor`, W, T, X, order)
}

c1_weight_1 <- function(W, T, X, inv, order) {
    .Call(`_VNM_c1_weight_1`, W, T, X, inv, order)
}

c1_weight_2 <- function(W, T, X, inv, order) {
    .Call(`_VNM_c1_weight_2`, W, T, X, inv, order)
}

c_weight_1 <- function(W, T, X, inv, dt, order) {
    .Call(`_VNM_c_weight_1`, W, T, X, inv, dt, order)
}

c_weight_2 <- function(W, T, X, inv, dt, order) {
    .Call(`_VNM_c_weight_2`, W, T, X, inv, dt, order)
}

D_weight_1 <- function(W, T, X, inv, order) {
    .Call(`_VNM_D_weight_1`, W, T, X, inv, order)
}

D_weight_2 <- function(W, T, X, inv, order) {
    .Call(`_VNM_D_weight_2`, W, T, X, inv, order)
}

M_weight_1 <- function(W, T, X, inv, dt, order, lambda) {
    .Call(`_VNM_M_weight_1`, W, T, X, inv, dt, order, lambda)
}

M_weight_2 <- function(W, T, X, inv, dt, order, lambda) {
    .Call(`_VNM_M_weight_2`, W, T, X, inv, dt, order, lambda)
}

rcpp_hello_world <- function() {
    .Call(`_VNM_rcpp_hello_world`)
}

