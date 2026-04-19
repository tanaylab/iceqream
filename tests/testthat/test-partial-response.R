test_that("escape_vars escapes regex metacharacters literally", {
    esc <- iceqream:::escape_vars("foo.bar")
    expect_false(grepl(esc, "foobar"))
    expect_true(grepl(esc, "foo.bar"))
    expect_false(grepl(esc, "fooxbar"))
})

test_that("escape_vars handles paren, question, pipe, plus, star, caret, dollar, brace", {
    for (s in c("a(b", "a?b", "a|b", "a+b", "a*b", "a^b", "a$b", "a{b", "a}b")) {
        esc <- iceqream:::escape_vars(s)
        expect_true(grepl(esc, s),
            info = sprintf("escape_vars(%s) = %s should match itself", s, esc)
        )
    }
})

