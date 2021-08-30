module Util
    export human_readable_size, parse_size, X, Y, Z, MIN, MAX

    # Aliases for common array indices
    const X = 1
    const Y = 2
    const Z = 3
    const MIN = 1
    const MAX = 2

    const ORDERS_OF_MAGNITUDE = "KMGTPEZY"

    """
    Convert `bytes` to a human-readable string.

    # Examples
    - 1024 → 1 KiB
    - 8589934592 → 8 GiB
    """
    function human_readable_size(bytes:: Integer) :: String
        order_of_magnitude = floor(Int, log(1024, bytes)) # Order of magnitude of bytes in base 1024 (0 for B, 1 for KiB, ...)
        order_of_magnitude = min(order_of_magnitude, length(ORDERS_OF_MAGNITUDE))
        bytes ÷= 1024^order_of_magnitude

        prefix = if order_of_magnitude == 0; ""
                 else ORDERS_OF_MAGNITUDE[order_of_magnitude] * "i" 
                 end

        return "$bytes $(prefix)B"
    end

    """
    Convert human-readable string representing a number of bytes to that number of bytes.

    # Examples
    - 1KiB → 1024
    - 8GiB → 8589934592
    """
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