process DOWNLOAD_GEO_ARCHIVE {
    label 'mbv_r'
    container r_image

    input:
    val r_image
    val timeout
    env 'GEO_ARCHIVE_RESOLVE_IP'

    output:
    path 'GSE333974_RAW.tar', emit: raw_archive

    script:
    """
    url='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE333nnn/GSE333974/suppl/GSE333974_RAW.tar'
    curl_args=(--fail --silent --show-error --location --connect-timeout 15 --max-time "$timeout")
    if [[ "\$GEO_ARCHIVE_RESOLVE_IP" != 'default' ]]; then
        curl_args+=(--resolve "ftp.ncbi.nlm.nih.gov:443:\$GEO_ARCHIVE_RESOLVE_IP")
    fi

    for attempt in 1 2 3; do
        echo "Downloading GSE333974_RAW.tar (attempt \$attempt/3)"
        if transfer=\$(curl "\${curl_args[@]}" \
            --output GSE333974_RAW.tar.part \
            --write-out '%{http_code} %{size_download} %{remote_ip}' "\$url"); then
            echo "GEO transfer (HTTP status, bytes, remote IP): \$transfer"
            if [[ "\$transfer" == '200 '* && -s GSE333974_RAW.tar.part ]]; then
                mv GSE333974_RAW.tar.part GSE333974_RAW.tar
                exit 0
            fi
        else
            curl_status=\$?
            if (( curl_status == 28 )); then
                echo "WARNING: GEO archive download attempt \$attempt/3 timed out (transfer limit: $timeout seconds; connection limit: 15 seconds). Consider increasing --timeout." >&2
            fi
            echo "GEO transfer failed (HTTP status, bytes, remote IP): \$transfer" >&2
        fi
        if (( attempt < 3 )); then
            sleep "\$((attempt * 2))"
        fi
    done

    echo 'Could not download the complete GSE333974_RAW.tar archive after 3 attempts.' >&2
    exit 1
    """
}
