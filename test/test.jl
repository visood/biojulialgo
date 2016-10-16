using Base.Test

function to_test_that(func::Function, message::ASCIIString)
    println("THAT $message is ", func() ? "TRUE" : "FALSE")
end

that(msg) = r::Test.Result ->
println("THAT $msg is", isa(r, Test.Success)? " TRUE" : " FALSE")

to_test = Test.with_handler

to_test_that(func::Function, msg::ASCIIString) =
    Test.with_handler( func, that(msg))
