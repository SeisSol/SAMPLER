module Util
    export human_readable_size

    const ORDERS_OF_MAGNITUDE = "KMGTPEZY"

    function human_readable_size(bytes:: Integer) :: String
        order_of_magnitude = floor(Int, log(1024, bytes)) # Order of magnitude of bytes in base 1024 (0 for B, 1 for KiB, ...)
        order_of_magnitude = min(order_of_magnitude, length(ORDERS_OF_MAGNITUDE))
        bytes ÷= 1024^order_of_magnitude

        prefix = if order_of_magnitude == 0; ""
                 else ORDERS_OF_MAGNITUDE[order_of_magnitude] * "i" 
                 end

        return "$bytes $(prefix)B"
    end

    function parse_size(size:: String) :: Integer
        if isnothing(match(r"^[0-9]+[KMGTPEZY]$", size))
            throw(ArgumentError("The given size string '$size' is invalid. Valid examples: 8G, 512M"))
        end

        magnitude_letter = size[end]
        magnitude = parse(Int, size[1:end-1])
        order_of_magnitude = findfirst(magnitude_letter, ORDERS_OF_MAGNITUDE)

        return magnitude * 1024^order_of_magnitude
    end
end